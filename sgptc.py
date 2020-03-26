from numpy import *
class particles(object):
    # def __new__(cls, fiefun = None, dt = 1e-10, num_ptc = 1, x = None, v= None, m = 1.674927211e-27, q = 1.602176634e-19):
    #     return super(particles, cls).__new__(cls)
    # Initialize fiefun is the field function which give E and B according to particle x and v.
    def __init__(self, fiefun = None, dt = 1e-10, num_ptc = 1, x = None, v= None, m = 1.674927211e-27, q = 1.602176634e-19, coordname = 'cylindrical'):
        self.dt = dt
        self.m = m
        self.q = q
        self.cd = coord(coordname)
        if x is None:
            x = array([[2.06, 0, 0]])
        if v is None:
            v = array([[9e3, 9e3, 2e3]])
        self.x = x
        self.v = v
        if fiefun == None:
            self.env = field((lambda x,v: zeros([2, x.shape[0], x.shape[1]])))
        else:
            self.env = field(fiefun)
        return

    #fields shape[2, m, 3],fields[0, :, :] is E, fields[1, :, :] is B. the second index is the label of particle num. the last one is three dimensions er, ez, ephi
    def solve_forces(self, x, v):
        self.env.evolution(x, v)#Since we don't solve fileds now,So this is useless now.
        self.fields = self.env.fields_return(x)#return fileds form x and v
        self.a = self.q / self.m * (self.fields[0, :, :] + cross(v, self.fields[1, :, :]))
        self.a += self.cd.dvdt_extra(self.x, self.v)
        # self.a[:, 0] = self.a[:, 0] + self.v[:, 2] * self.v[:, 2] / self.x[:, 0]#centrifugal force
        # self.a[:, 2] = self.a[:, 2] - self.v[:, 0] * self.v[:, 2] / self.x[:, 0]        
        return self.a
        
    def evo_v(self):
        self.v += self.a * self.dt
        return

    #Since in boris method, we push particle with E and B individually,So we need find a just from E
    def solve_forces_e(self, x, v):
        self.env.evolution(x, v)
        self.fields = self.env.fields_return(x)
        self.ae = self.q / self.m * (self.fields[0, :, :])
        self.ae += self.cd.dvdt_extra(self.x, self.v)
        # self.ae[:, 0] = self.ae[:, 0] + self.v[:, 2] * self.v[:, 2] / self.x[:, 0]
        # self.ae[:, 2] = self.ae[:, 2] - self.v[:, 0] * self.v[:, 2] / self.x[:, 0]

    # push half dt
    def evo_v_boris_e(self):
        self.v += self.ae * self.dt / 2
        
    
    def evo_v_boris_b(self):
        t = self.q / self.m * self.fields[1, :, :] * self.dt / 2
        vp = self.v + cross(self.v, t)
        s = 2 * t / (1 + (t * t).sum())
        self.v += cross(vp, s)
        
    
    def evo_x(self):
        v = self.v.copy()
        self.x += self.cd.dvdt_extra(self.x, v)
        self.x += self.dt * v
        return

    def boris(self):
        self.solve_forces_e(self.x, self.v)
        self.evo_v_boris_e()
        self.evo_v_boris_b()
        self.evo_v_boris_e()
        self.evo_x()
        return
        
    def leap_frog(self):
        self.solve_forces(self.x, self.v)
        self.evo_v()
        self.evo_x()
        return
        
    def rk4(self):
        t4 = array([self.dt/2, self.dt/2, self.dt, self.dt])
        mul = [1/6.0, 1/3.0, 1/3.0, 1/6.0]
        # fieldstemp = self.env.fields_return(self.x)
        # atemp = [self.q * (fieldstemp[:, :, 0] + cross(self.v, fieldstemp[:, :, 1])) / self.m]
        xtemp = [self.x]        
        vtemp = [self.v]
        vtempd = self.v.copy()
        vtempd[:, 2] /= xtemp[0][:, 0]
        vd = vtempd.copy() * mul[0]
        atemp = []
        for i in range(3):
            self.solve_forces(xtemp[i], vtemp[i])
            atemp.append(self.a)
            vtemp.append(atemp[i] * t4[i] + self.v)
            xtemp.append(vtempd * t4[i] + self.x)
            vtempd = vtemp[i].copy()
            vtempd[:, 2] /= xtemp[i][:, 0] 
            vd += vtempd * mul[i+1]
        atemp.append(self.solve_forces(xtemp[3], vtemp[3]))
        a = atemp[0] / 6 + atemp[1] / 3 + atemp[2] / 3 + atemp[3] / 6
        self.v += a * self.dt
        self.x += vd * self.dt 
        return

    # def __leap_frog_initial__(self):
    #     self.leap_frog_initial.__func__.__code__ = (lambda x:None).__code__

    # def __get_fields_function__(self, fields_return = (lambda x: rollaxis(array([zeros((x.shape)), ones((x.shape))]), 0,3))):
    #     self.field_return = fields_return
    #     return

    
class field(object):
    #initialize,If the fields_function was given then use it.If not use the default one:E=0 everywhere Bphi=1 B_others = 0 everywhere
    def __init__(self,fields_return = (lambda x: array([zeros((x.shape)), rollaxis(array([zeros(x.shape[0]), zeros(x.shape[0]), ones(x.shape[0])]), 0, 2)]))):
        self.fields_return = fields_return
        return
    def evolution(self, x, v):
        return

class coord(object):
    def __init__(self, coord = 'cylindrical'):
        self.coordname = coord
        coord += ';'
        cmd = 'self.dvdt_extra = self.dvdt_extra_' + coord + 'self.dxdt_extra_ = self.dxdt_extra_' + coord
        # self.vr_dvdt = x.copy()
        # self.vr_drdt = x.copy()
        exec(cmd)
        
    def dvdt_extra_cylindrical(self, x, v):
        aextra = zeros(x.shape)
        aextra[:, 0] = v[:, 2] * v[:, 2] / x[:, 0]#centrifugal force
        aextra[:, 2] = - v[:, 0] * v[:, 2] / x[:, 0]        
        return aextra
    def dxdt_extra_cylindrical(self, x, v):
        vextra = zeros(v.shape)
        vextra[:, 2] = v[:, 2] / x[:, 0] -v[:, 2]
        
        return 0

    def dvdt_extra_cartesian(self, x, v):
        
        return 0

    def dxdt_extra_cartesian(self, ):
        
        return 0

    def jacobian_cylindrical():
        
        return

    def jacobian_cartesian():

        return

def cross(v, b):
    re = array((v[:, 1] * b[:, 2] - v[:, 2] * b[:, 1], v[:, 2] * b[:, 0] - v[:, 0] * b[:, 2], v[:, 0] * b[:, 1] - v[:, 1] * b[:, 0])).T
    return re
