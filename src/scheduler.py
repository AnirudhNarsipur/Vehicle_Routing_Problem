from __future__ import annotations

from dataclasses import dataclass

from docplex.cp.model import *


@dataclass
class Config:
    n_weeks: int
    n_days: int
    n_days_in_week: int  # computed as n_days // n_weeks
    n_employees: int
    n_shifts: int
    n_intervals_in_day: int
    min_shifts: list[list[int]]  # demand per day per shift
    min_daily: int
    employee_min_daily: int  # except on off-shifts
    employee_max_daily: int
    employee_min_weekly: int
    employee_max_weekly: int
    employee_max_consecutive_night_shifts: int
    employee_max_total_night_shifts: int

    @staticmethod
    def load(f) -> Config:
        """
        Reads in a file of the following format:
        Business_numWeeks: 1
        Business_numDays: 7
        Business_numEmployees: 14
        Business_numShifts: 4
        Business_numIntervalsInDay: 24
        Business_minDemandDayShift: 0 1 3 2 0 1 3 2 0 1 3 2 0 1 5 2 0 1 5 2 0 1 5 2 0 1 5 2
        Business_minDailyOperation: 60
        Employee_minConsecutiveWork: 4
        Employee_maxDailyWork: 8
        Employee_minWeeklyWork: 20
        Employee_maxWeeklyWork: 40
        Employee_maxConsecutiveNigthShift: 1
        Employee_maxTotalNigthShift: 2
        :param f: path to the file
        :return: a dictionary with the parameters
        """
        with open(f, "r") as fl:
            lines = fl.readlines()
            params = {}
            for line in lines:
                line = line.strip()
                if line.startswith("#"):
                    continue
                if line.startswith("Business_"):
                    key, value = line.split(":")
                    if key != "Business_minDemandDayShift":
                        params[key] = int(value)
                    else:
                        params[key] = [int(x) for x in value.split()]
                elif line.startswith("Employee_"):
                    key, value = line.split(":")
                    params[key] = int(value)
                else:
                    raise ValueError("Invalid line in file")
        n_weeks = params["Business_numWeeks"]
        n_days = params["Business_numDays"]
        n_days_in_week = n_days // n_weeks
        n_employees = params["Business_numEmployees"]
        n_shifts = params["Business_numShifts"]
        n_intervals_in_day = params["Business_numIntervalsInDay"]
        min_shifts = []
        for i in range(0, n_days * n_shifts, n_shifts):
            min_shifts.append(params["Business_minDemandDayShift"][i : i + n_shifts])
        min_daily = params["Business_minDailyOperation"]
        employee_min_daily = params["Employee_minConsecutiveWork"]
        employee_max_daily = params["Employee_maxDailyWork"]
        employee_min_weekly = params["Employee_minWeeklyWork"]
        employee_max_weekly = params["Employee_maxWeeklyWork"]
        employee_max_consecutive_night_shifts = params[
            "Employee_maxConsecutiveNigthShift"
        ]
        employee_max_total_night_shifts = params["Employee_maxTotalNigthShift"]
        return Config(
            n_weeks,
            n_days,
            n_days_in_week,
            n_employees,
            n_shifts,
            n_intervals_in_day,
            min_shifts,
            min_daily,
            employee_min_daily,
            employee_max_daily,
            employee_min_weekly,
            employee_max_weekly,
            employee_max_consecutive_night_shifts,
            employee_max_total_night_shifts,
        )


Schedule = list[list[tuple[int, int]]]


@dataclass
class Solution:
    is_solution: bool
    n_fails: int
    schedule: Schedule


class Scheduler:
    def __init__(self, config: Config):
        self.config = config
        self.model = CpoModel()
        self.build_constraints()

    def build_constraints(self):
        self.shifts = []
        self.hours = []
        # domains cover min and max daily works
        # shift length is covered by the shift lengths for 24 hour, 3 shift days
        # this is not true in general, but we leave it out for optimization
        for _ in range(self.config.n_employees):
            self.shifts.append(
                integer_var_list(
                    self.config.n_days, min=0, max=self.config.n_shifts - 1
                )
            )
            self.hours.append(
                integer_var_list(
                    self.config.n_days,
                    domain=[
                        0,
                        (
                            self.config.employee_min_daily,
                            self.config.employee_max_daily,
                        ),
                    ],
                )
            )

        # non-off shifts must not be 0 hours, and off shifts must be 0 hours
        for e in range(self.config.n_employees):
            for n in range(self.config.n_days):
                self.model.add(if_then(self.shifts[e][n] == 0, self.hours[e][n] == 0))
                self.model.add(if_then(self.shifts[e][n] != 0, self.hours[e][n] > 0))

        # business needs
        for n in range(self.config.n_days):
            self.model.add(
                sum([self.hours[e][n] for e in range(self.config.n_employees)])
                >= self.config.min_daily
            )
            for shift in range(self.config.n_shifts):
                self.model.add(
                    count(
                        [self.shifts[e][n] for e in range(self.config.n_employees)],
                        shift,
                    )
                    >= self.config.min_shifts[n][shift]
                )

        # training requirements
        for e in range(self.config.n_employees):
            self.model.add(all_diff(self.shifts[e][: self.config.n_shifts]))

        # weekly work
        for e in range(self.config.n_employees):
            for w in range(self.config.n_weeks):
                work_in_week = sum(
                    [
                        self.hours[e][d]
                        for d in range(
                            w * self.config.n_days_in_week,
                            min(
                                (w + 1) * self.config.n_days_in_week, self.config.n_days
                            ),
                        )
                    ]
                )
                self.model.add(work_in_week >= self.config.employee_min_weekly)
                self.model.add(work_in_week <= self.config.employee_max_weekly)

        # TODO: night shift rules

    def solve(self) -> Solution:
        print(self.config)
        solution = self.model.solve()
        n_fails = solution.get_solver_info(CpoSolverInfos.NUMBER_OF_FAILS)
        if not solution:
            return Solution(false, n_fails, None)
        schedule = self.construct_schedule(solution)
        return Solution(true, n_fails, schedule)

    def construct_schedule(self, solution: CpoSolveResult) -> Schedule:
        shift_starts = {}
        for i in range(1, self.config.n_shifts):
            shift_starts[i] = (
                (i - 1) * self.config.n_intervals_in_day // (self.config.n_shifts - 1)
            )
        schedule = []
        for i in range(self.config.n_employees):
            employee_schedule = []
            employee_shifts = self.shifts[i]
            employee_hours = self.hours[i]
            for j in range(self.config.n_days):
                shift = solution[employee_shifts[j]]
                hours = solution[employee_hours[j]]
                if shift == 0:  # off shift
                    employee_schedule.append((-1, -1))
                else:
                    employee_schedule.append(
                        (shift_starts[shift], shift_starts[shift] + hours)
                    )
            schedule.append(employee_schedule)

        return schedule

    @staticmethod
    def from_file(f) -> Scheduler:
        config = Config.load(f)
        return Scheduler(config)
