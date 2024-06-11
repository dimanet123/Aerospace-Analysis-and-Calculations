import numpy as np

class Planet:
    def __init__(self, mu, R, eps, We, sol_activity, meters = True, sol_activity_bool = True):
        if meters == True:
            self.mu = mu
            self.R = R
            self.eps = eps
            self.We = We
            self.sol_activity = sol_activity
            self.sol_activity_bool = sol_activity_bool
        else:
            self.mu = mu*1e9
            self.R = R*1e3
            self.eps = eps*1e15
            self.We = We
            self.sol_activity = sol_activity
            self.sol_activity_bool = sol_activity_bool
        
    def calc_atm(self, h):
        if self.sol_activity_bool == True:
            """ Расчёт плотности атмосферы на высоте h с учётом солнечной активности s. """
            c = [-19.6120560797504, -4.37607377008974E-02, 1.03643542927751E-05, 3.19927015574673E-08,
                 -1.84888682338092E-06, 1.08392130394651E-04, 1.64994623915598E-07,
                 -4.50548229523859E-10, 6.97982431033395E-10, -2.00112208505143E-07,
                 -2.85143658777182E-10, 9.3931325534369E-13, -1.99285975233419E-12]
            h = h/ 1000
            y1 = h - 120
            y2 = 0
            if h < 225:
                y2 = h - 225

            l = c[0] + c[1] * y1 + c[2] * y1**2 + c[3] * y1**3 + c[4] * y2**3
            l += self.sol_activity * (c[5] * y1 + c[6] * y1**2 + c[7] * y1**3 + c[8] * y2**3)
            l += self.sol_activity**2 * (c[9] * y1 + c[10] * y1**2 + c[11] * y1**3 + c[12] * y2**3)

            return np.exp(l)
        else:
            a = np.array([-6.3759, -7.3012, -1.1817])
            b = np.array([-0.4754, -0.0096, -0.0068, -0.0120, 0.0042])
            c = np.array([0.1803, 0.0872, -0.0153, 0.0145, 0])

            # Функция вычисления плотности воздуха
            # высота h задаётся в километрах
            # результат -- плотность кг/м3
            h = h / 1000
            x = h/50.0-1
            sa = a[0] + a[1]*x + a[2]*x*x
            sbc = np.sum((b[i]*np.cos((i+1)*np.pi*x) + c[i] *
                          np.sin((i+1)*np.pi*x) for i in range(5)))
            return np.exp(sa + sbc)
    