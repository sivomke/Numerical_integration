from blackHole import TargetFunction, middle_rect
import scipy
import scipy.integrate


def third_derivative(a, b, x):
    return scipy.exp(-a*x)*((3*pow(a,2)*b-pow(a,3))*scipy.sin(b*x)+(3*pow(a,2)*b-pow(b,3))*scipy.cos(b*x))


# returns a tuple - step h and auxiliary variable i
def find_step(a, b, left, right, k):
    w = TargetFunction(a, b)
    extremums = []
    extremums.extend([left, right])
    i = k
    x = scipy.arctan((pow(w.b,3)-3*pow(w.a,2)*w.b)/(3*pow(w.a,2)*w.b-pow(w.a, 3)))/w.b + i*scipy.pi/w.b
    extremums.append(x)
    while True:
        i += 1
        x = scipy.arctan((pow(w.b,3)-3*pow(w.a,2)*w.b)/(3*pow(w.a,2)*w.b-pow(w.a,3)))/w.b + i*scipy.pi/w.b
        if x > right:
            break
        else:
            extremums.append(x)
    values = list(map(lambda x: abs(third_derivative(w.a, w.b, x)), extremums))
    m = max(values)  # max of third derivative of target func on given interval
    h = right - left
    eps = m*pow(h, 2)*(right - left)/24
    while eps > 0.5*pow(10, -5):
        h /= 2
        eps = m * pow(h, 2) * (right - left) / 24
    return h, i


def composite_rect(step, left, right, f):
    count = int((right - left)/step)
    tmr_left = left
    s = 0
    for i in range(0, count):
        s += middle_rect(tmr_left, tmr_left + step, f)
        tmr_left += step
    return s


def main():
    g = TargetFunction(2, 5)
    left = 0
    right = scipy.pi / g.b
    aux = find_step(g.a, g.b, left, right, 0)
    h = aux[0]
    print("Step of currant quadrature: ", '{:.5f}'.format(h))
    k = aux[1] # auxiliary variable to reduce num of calculation
    s = composite_rect(h, left, right, g.value)
    while True:
        left += scipy.pi / g.b
        right += scipy.pi / g.b
        aux = find_step(g.a, g.b, left, right, k)
        h = aux[0]
        print("Step of currant quadrature: ", '{:.5f}'.format(h))
        k = aux[1]
        quad = composite_rect(h, left, right, g.value)
        if abs(quad) <= g.tolerance / 2:
            break
        else:
            s += quad
    print("\n")
    print("Approximate integral value: ", format(s, '.5f'))
    print("Tolerance: ", g.tolerance)

# main()






