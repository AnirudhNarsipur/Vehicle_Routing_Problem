import os
from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

from docplex.cp.config import context

from src.scheduler import Schedule, Scheduler

context.solver.agent = "local"
context.solver.local.execfile = os.environ["CP_SOLVER_EXEC"]

parser = ArgumentParser()
parser.add_argument("input_file", type=str)

# lifted right out of the Java support code. Janky af
def visualize(schedule: Schedule):
    for e in range(len(schedule)):
        print("E" + str(e + 1) + ": ", end="")
        if e < 9:
            print(" ", end="")
        for d in range(len(schedule[0])):
            for i in range(24):
                if i % 8 == 0:
                    print("|", end="")
                if (
                    schedule[e][d][0] != schedule[e][d][1]
                    and i >= schedule[e][d][0]
                    and i < schedule[e][d][1]
                ):
                    print("+", end="")
                else:
                    print(".", end="")

            print("|", end="")
        print(" ")


def main(args):
    input_file = Path(args.input_file)
    filename = input_file.name
    scheduler = Scheduler.from_file(args.input_file)
    start_time = datetime.now()
    solution = scheduler.solve()
    serialized_schedule = []
    for employee_schedule in solution.schedule:
        for start, end in employee_schedule:
            serialized_schedule.append(str(start))
            serialized_schedule.append(str(end))
    serialized_schedule = " ".join(serialized_schedule)
    end_time = datetime.now()
    delta = round((end_time - start_time).total_seconds() * 100) / 100
    visualize(solution.schedule)
    print(
        f'{{"Instance": "{filename}", "Time": {delta}, "Result": {solution.n_fails}, "Solution": "{serialized_schedule}"}}'
    )


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)