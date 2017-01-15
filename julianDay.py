# Handy julian date converter! Input list in the form
# [year, month, day,hour, minute, second]


def toJulian(time):
    a = int(time[0] / 100)
    b = int(a / 4)
    c = 2 - a + b
    e = int(365.25 * (time[0] + 4716))
    f = int(30.6001 * (time[1] + 1))
    preJulian = c + time[2] + e + f - 1524.5
    hours = (time[3] + time[4] / float(60) + time[5] / float(3600)) / float(24)
    julian = preJulian + hours
    return julian


print toJulian([2015, 7, 26, 4, 19, 54])
print toJulian([2015, 7, 26, 4, 36, 51])
print toJulian([2015, 7, 26, 4, 53, 25])
