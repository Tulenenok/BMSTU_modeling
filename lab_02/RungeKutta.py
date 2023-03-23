class RungeKutta:
    def __init__(self, dIdt, dUdt, t0, I0, U0, tm, h):
        self.dIdt = dIdt
        self.dUdt = dUdt
        self.t0 = t0
        self.I0 = I0
        self.U0 = U0
        self.tm = tm
        self.h = h

    def order_accuracy_1(self):
        last_i = (self.t0, self.I0)
        last_u = (self.t0, self.U0)

        i_res = []
        u_res = []

        t = self.t0
        while t <= self.tm:
            I = last_i[1] + self.h * self.dIdt(last_i[0], last_u[1], last_i[1])
            U = last_u[1] + self.h * self.dUdt(last_i[0], last_u[1], last_i[1])

            i_res.append(last_i)
            u_res.append(last_u)

            last_i = (t, I)
            last_u = (t, U)

            t += self.h

        return i_res, u_res

    def order_accuracy_2(self):
        beta = 1
        last_i = (self.t0, self.I0)
        last_u = (self.t0, self.U0)

        i_res = []
        u_res = []

        t = self.t0
        while t <= self.tm:
            I1 = self.h * self.dIdt(last_i[0], last_u[1], last_i[1])
            U1 = self.h * self.dUdt(last_i[0], last_u[1], last_i[1])

            I2 = self.h * self.dIdt(last_i[0] + self.h / (2 * beta), last_u[1] + U1 / (2 * beta), last_i[1] + I1 / (2 * beta))
            U2 = self.h * self.dUdt(last_i[0] + self.h / (2 * beta), last_u[1] + U1 / (2 * beta), last_i[1] + I1 / (2 * beta))

            I = last_i[1] + ((1 - beta) * I1 + beta * I2)
            U = last_u[1] + ((1 - beta) * U1 + beta * U2)

            i_res.append(last_i)
            u_res.append(last_u)

            last_i = (t, I)
            last_u = (t, U)

            t += self.h

        return i_res, u_res

    def order_accuracy_4(self):
        beta = 1
        last_i = (self.t0, self.I0)
        last_u = (self.t0, self.U0)

        i_res = []
        u_res = []

        t = self.t0
        while t <= self.tm:
            I1 = self.h * self.dIdt(last_i[0], last_u[1], last_i[1])
            U1 = self.h * self.dUdt(last_i[0], last_u[1], last_i[1])

            I2 = self.h * self.dIdt(last_i[0] + self.h / 2, last_u[1] + U1 / 2, last_i[1] + I1 / 2)
            U2 = self.h * self.dUdt(last_i[0] + self.h / 2, last_u[1] + U1 / 2, last_i[1] + I1 / 2)

            I3 = self.h * self.dIdt(last_i[0] + self.h / 2, last_u[1] + U2 / 2, last_i[1] + I2 / 2)
            U3 = self.h * self.dUdt(last_i[0] + self.h / 2, last_u[1] + U2 / 2, last_i[1] + I2 / 2)

            I4 = self.h * self.dIdt(last_i[0] + self.h, last_u[1] + U3, last_i[1] + I3)
            U4 = self.h * self.dUdt(last_i[0] + self.h, last_u[1] + U3, last_i[1] + I3)

            I = last_i[1] + (I1 + 2 * I2 + 2 * I3 + I4) / 6
            U = last_u[1] + (U1 + 2 * U2 + 2 * U3 + U4) / 6

            i_res.append(last_i)
            u_res.append(last_u)

            last_i = (t, I)
            last_u = (t, U)

            t += self.h

        return i_res, u_res

    def __call__(self, order_accuracy):
        if order_accuracy == 1:
            return self.order_accuracy_1()
        elif order_accuracy == 2:
            return self.order_accuracy_2()
        elif order_accuracy == 4:
            return self.order_accuracy_4()
        else:
            print('Неверно введен порядок точности')
            return None, None
