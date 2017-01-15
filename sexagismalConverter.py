# Another handy converter!


def toDecimal(sexagismal):
    if sexagismal[0] > 1:
        time = sexagismal[0]
        time += (sexagismal[1] / float(60))
        time += (sexagismal[2] / float(3600))
    else:
        time = sexagismal[0]
        time -= (sexagismal[1] / float(60))
        time -= (sexagismal[2] / float(3600))
    return time

# Example use
print toDecimal([-7, 20, 24.3])
print toDecimal([17, 6, 26.13])
