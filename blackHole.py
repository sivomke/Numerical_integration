import scipy
import scipy.integrate


# defining target function

class TargetFunction:
    def __init__(self, a, b):
        self.a = a
        self.b = b
        self.tolerance = pow(10, -5)

    def value(self, x):
        return scipy.exp(-self.a*x)*scipy.sin(self.b*x)


# implements middle rectangle method
def middle_rect(left, right, f):
    return (right-left)*f(0.5*(right+left))


def adaptive(left, right, f, tolerance):
    s = 0
    tmr_left = left
    h = (right - left)/10
    steps = []
    while True:
        while True:
            rect = middle_rect(tmr_left, tmr_left + h, f)
            rect_comp = middle_rect(tmr_left, tmr_left + 0.5*h, f) + middle_rect(tmr_left + 0.5*h, tmr_left + h, f)
            eps = (rect_comp - rect) / 3
            if abs(eps) <= h*tolerance/2*(right-left):
                break
            else:
                h *= 0.5
        steps.append(h)
        tmr_left += h
        s += rect_comp
        h *= 2
        if tmr_left == right:
            break
        if (tmr_left + h) > right:
            h = right - tmr_left
    print("adaptive steps: ")
    print(', '.join('{:.5f}'.format(el) for el in steps))
    return s


def main():
    g = TargetFunction(2, 5)
    left = 0
    right = scipy.pi/g.b
    print()
    s = adaptive(left, right, g.value, g.tolerance)
    while True:
        left += scipy.pi/g.b
        right += scipy.pi/g.b
        quad = adaptive(left, right, g.value, g.tolerance)
        if abs(quad) <= g.tolerance/2:
            break
        else:
            print("")
            s += quad
    print("\n")
    print("Approximate integral value: ", format(s, '.5f'))
    print("Tolerance: ", g.tolerance)


main()



