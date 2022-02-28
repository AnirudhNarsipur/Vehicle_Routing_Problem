# Function read in a file of the following format:
def parse(fl_path : str) -> dict:
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
    :param fl_path: path to the file
    :return: a dictionary with the parameters
    """
    with open(fl_path, 'r') as fl:
        lines = fl.readlines()
        params = {}
        for line in lines:
            line = line.strip()
            if line.startswith('#'):
                continue
            if line.startswith('Business_'):
                key, value = line.split(':')
                if key != "Business_minDemandDayShift":
                    params[key] = int(value)
                else:
                    params[key] = [int(x) for x in value.split()]
            elif line.startswith('Employee_'):
                key, value = line.split(':')
                params[key] = int(value)
            else:
                raise ValueError('Invalid line in file')
    return params