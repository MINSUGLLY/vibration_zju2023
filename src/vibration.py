import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from scipy import signal
from scipy import interpolate
import matplotlib.animation as ani
import webbrowser

# 边界条件编号  0自由 1简支 2固支  

class Beam:
    def __init__(self, name = 'beam', l = 1, b = lambda x : 0.1, h = lambda x : 0.1, rho = lambda x : 7850, E = lambda x : 2e11, f = lambda x, t : 10000 * np.sin(2 * np.pi * t), bdL = 2, bdR = 0):
        assert bdL + bdR >= 2
        self.l = l
        self.name = name
        self.b = b
        self.h = h
        self.rho = rho
        self.E = E
        self.f = f
        self.bdL = bdL
        self.bdR = bdR

        self.S = lambda x : self.b(x) * self.h(x)
        self.I = lambda x : self.b(x) * self.h(x) ** 3 / 12
    
class FEM(Beam):
    def __init__(self,name = 'beam' , l=1, b=lambda x: 0.1,h=lambda x: 0.1, rho=lambda x: 7850, E=lambda x: 2e11, f = lambda x, t : 10000 * np.sin(2 * np.pi * t), bdL=2, bdR=0, x_step=0.01, t_end=0.1, t_step=1e-5):
        super().__init__(name, l, b, h, rho, E, f, bdL, bdR)
        self.x_step = x_step
        self.xs = np.arange(0, self.l+x_step, x_step)
        self.t_end = t_end
        self.t_step = t_step
        self.ts = np.arange(0, t_end+t_step, t_step)
        
        self.N = [lambda x : 1 - 3*x**2/x_step**2 + 2*x**3/x_step**3, lambda x : x - 2*x**2/x_step + x**3/x_step**2, lambda x : 3*x**2/x_step**2 - 2*x**3/x_step**3, lambda x : -x**2/x_step + x**3/x_step**2]
        self.N__ = [lambda x : -6/x_step**2 + 12*x/x_step**3, lambda x : -4/x_step + 6*x/x_step**2, lambda x : 6/x_step**2 - 12*x/x_step**3, lambda x : -2/x_step + 6*x/x_step**2]
        self.elementNum = int(self.l/x_step) if abs(self.l/x_step - int(self.l/x_step)) < 1e-9 else int(self.l/x_step) + 1
        
        M = np.zeros((4*self.elementNum, 4*self.elementNum))
        K = np.zeros((4*self.elementNum, 4*self.elementNum))
        self.vL = np.zeros(self.elementNum)
        for i in range(self.elementNum):
            self.vL[i] = i * self.x_step
            for j in range(4):
                for k in range(4):
                    M[4*i+j, 4*i+k] = integrate.quad(lambda x : self.rho(self.vL[i]+x)*self.S(self.vL[i]+x)*self.N[j](x)*self.N[k](x), 0, x_step)[0]
                    K[4*i+j, 4*i+k] = integrate.quad(lambda x : self.E(self.vL[i]+x)*self.I(self.vL[i]+x)*self.N__[j](x)*self.N__[k](x), 0, x_step)[0]
                
        self.DoF = 2*(self.elementNum+1) - self.bdL - self.bdR
        self.beta = np.zeros((4*self.elementNum, self.DoF))
        if self.bdL == 0:
            self.beta[0][0] = 1
            self.beta[1][1] = 1
        elif self.bdL ==1:
            self.beta[1][0] = 1
        
        if self.bdR == 0:
            self.beta[-1][-1] = 1
            self.beta[-2][-2] = 1
        elif self.bdR ==1:
            self.beta[-1][-1] = 1
        
        for i in range(2, 4*self.elementNum-2):
            if (i-2)%4 == 0 or (i-2)%4 == 2:
                self.beta[i][2*int((i-2)/4)+2-bdL] = 1
            else:
                self.beta[i][2*int((i-2)/4)+3-bdL] = 1
        F = lambda t, i : integrate.quad(lambda x : self.f(self.vL[i]+x, t) * self.N[i](x), 0, self.x_step)[0]
        self.M = np.array(np.mat(self.beta).transpose() * np.mat(M) * np.mat(self.beta))
        self.K = np.array(np.mat(self.beta).transpose() * np.mat(K) * np.mat(self.beta))
        self.Q = lambda t, i : float(np.mat(self.beta[:, i]) * np.mat([F(t, j%4) for j in range(4*self.elementNum)]).transpose())
        
        D = np.array(np.linalg.inv(np.mat(self.K)) * np.mat(self.M))
        eig_value, self.phi = np.linalg.eig(D)
        self.omega = 1/eig_value**0.5
        print('系统固有频率为：')
        print(self.omega)
        
    
    
    def ModeSuperpositionMethod(self):        
        M_p = []
        K_p = []
        Q_p = lambda t, i : float(np.mat(self.phi[:, i]) * np.mat([self.Q(t, j) for j in range(self.DoF)]).transpose())
        for i in range(self.DoF):
            phi_i = np.mat(self.phi[:, i]).transpose()
            M_p.append(float(np.array(phi_i.transpose() * np.mat(self.M) * phi_i)))
            K_p.append(float(np.array(phi_i.transpose() * np.mat(self.K) * phi_i)))
        M_p = np.array(M_p)
        K_p = np.array(K_p)
        
        h_p = lambda t, i : np.array(1/(M_p[i]*self.omega[i]) * np.sin(self.omega[i]*t))
        
        h_p_discrete = []
        Q_p_discrete = []
        for i in range(self.DoF):
            h_p_discrete.append(h_p(self.ts, i))
            Q_p_discrete_i = []
            for t in self.ts:
                Q_p_discrete_i.append(Q_p(t, i))
            Q_p_discrete.append(Q_p_discrete_i)
        h_p_discrete = np.array(h_p_discrete)
        Q_p_discrete = np.array(Q_p_discrete)

        q_p = []
        for i in range(self.DoF):
            q_p.append(signal.convolve(h_p_discrete[i], Q_p_discrete[i])[0:len(self.ts)]*self.t_step) 
        q_p = np.array(q_p)
        q = np.array(np.mat(self.phi) * np.mat(q_p))
        
        w_theta = np.array(np.mat(self.beta) * np.mat(q))
        
        deleteList = []
        for i in range(self.elementNum-1):
            deleteList.append(4*i+4)
            deleteList.append(4*i+5)
        w_theta = np.delete(w_theta, deleteList, axis=0)
        
        interpolate_xs = np.zeros(2*(self.elementNum+1))
        for i in range(self.elementNum+1):
            interpolate_xs[2*i] = self.xs[i]
            interpolate_xs[2*i+1] = self.xs[i]

        plot_xs = np.arange(0, self.l + 1e-8, 0.01)
        plot_ys = []
        for i in range(len(self.ts)):
            interpolate_ys = w_theta[:, i]
            interpolant = interpolate.KroghInterpolator(interpolate_xs, interpolate_ys)
            plot_ys.append(interpolant(plot_xs))
        
        return plot_xs, plot_ys
    
    
        
    def Newmark_beta(self, beta = 0.25, gamma = 0.5):
        a0 = 1/(beta*self.t_step**2)
        a1 = gamma/(beta*self.t_step)
        a2 = 1/(beta*self.t_step)
        a3 = 1/(2*beta) - 1
        a4 = gamma/beta - 1
        a5 = self.t_step/2*(gamma/beta - 1)
        a6 = self.t_step*(1-gamma)
        a7 = gamma*self.t_step
        
        q = np.zeros((self.DoF, len(self.ts)))
        dq = np.zeros((self.DoF, 1))
        ddq = np.array(np.linalg.inv(self.M) * np.mat(np.array([self.Q(0, i) for i in range(self.DoF)])).transpose()).reshape(self.DoF, 1)
        K_ = a0*self.M + self.K
        
        
        for t in range(1, len(self.ts)):
            h = a0*q[:, t-1].reshape((self.DoF, 1)) + a2*dq +a3*ddq
            f_ = lambda t, i : self.Q(t, i) + np.dot(self.M[i, :], h)
            q[:, t] = np.array(np.mat(np.linalg.inv(K_)) * np.mat([f_(t*self.t_step, i) for i in range(self.DoF)])).reshape(self.DoF,)
            ddq_old = ddq
            ddq = a0*(q[:, t]-q[:, t-1]).reshape((self.DoF, 1)) - a2*dq -a3*ddq
            dq = dq + a6*ddq_old + a7*ddq
        
        w_theta = np.array(np.mat(self.beta) * np.mat(q))
        
        deleteList = []
        for i in range(self.elementNum-1):
            deleteList.append(4*i+4)
            deleteList.append(4*i+5)
        w_theta = np.delete(w_theta, deleteList, axis=0)
        
        interpolate_xs = np.zeros(2*(self.elementNum+1))
        for i in range(self.elementNum+1):
            interpolate_xs[2*i] = self.xs[i]
            interpolate_xs[2*i+1] = self.xs[i]

        plot_xs = np.arange(0, self.l + 1e-8, 0.01)
        plot_ys = []
        for i in range(len(self.ts)):
            interpolate_ys = w_theta[:, i]
            interpolant = interpolate.KroghInterpolator(interpolate_xs, interpolate_ys)
            plot_ys.append(interpolant(plot_xs))
        
        return plot_xs, plot_ys
    
    def Wilson_theta(self, theta = 1.5):
        q = np.zeros((self.DoF, len(self.ts)))
        dq = np.zeros((self.DoF, 1))
        ddq = np.array(np.linalg.inv(self.M) * np.mat(np.array([self.Q(0, i) for i in range(self.DoF)])).transpose()).reshape(self.DoF, 1)
        
        K_ = self.K + 6*self.M/(theta*self.t_step)**2
        
        for t in range(1, len(self.ts)):
            p_ = np.array([theta*self.Q(t*self.t_step, i) + (1-theta)*self.Q((t-1)*self.t_step, i) for i in range(self.DoF)]).reshape((self.DoF, 1)) + np.array(np.mat(self.M)*(6*(q[:, t-1]).reshape((self.DoF, 1))/(theta*self.t_step)**2 + 6*dq/(theta*self.t_step) + 2*ddq)).reshape((self.DoF, 1))
            delta_q = np.array(np.mat(np.linalg.inv(K_)) * np.mat(p_)).reshape((self.DoF, 1))
            
            ddq_old = ddq
            ddq = 6/(theta**3*self.t_step**2)*(delta_q - q[:, t-1].reshape((self.DoF, 1))) - 6/(theta**2*self.t_step)*dq + (1 - 3/theta)*ddq
            dq_old = dq
            dq = dq + self.t_step/2*(ddq+ddq_old)
            q[:, t] = q[:, t-1] + (self.t_step*dq_old + (self.t_step**2)/6*(ddq+2*ddq_old)).reshape(self.DoF,)
        
        w_theta = np.array(np.mat(self.beta) * np.mat(q))
        
        deleteList = []
        for i in range(self.elementNum-1):
            deleteList.append(4*i+4)
            deleteList.append(4*i+5)
        w_theta = np.delete(w_theta, deleteList, axis=0)
        
        interpolate_xs = np.zeros(2*(self.elementNum+1))
        for i in range(self.elementNum+1):
            interpolate_xs[2*i] = self.xs[i]
            interpolate_xs[2*i+1] = self.xs[i]

        plot_xs = np.arange(0, self.l + 1e-8, 0.01)
        plot_ys = []
        for i in range(len(self.ts)):
            interpolate_ys = w_theta[:, i]
            interpolant = interpolate.KroghInterpolator(interpolate_xs, interpolate_ys)
            plot_ys.append(interpolant(plot_xs))
        
        return plot_xs, plot_ys
            
    
    def draw(self):
        x1, y1 = self.ModeSuperpositionMethod()
        y1 = np.array(y1)[np.arange(0, len(y1), 10), :]
        print('模态叠加法 完成')
        x2, y2 = self.Newmark_beta()
        y2 = np.array(y2)[np.arange(0, len(y2), 10), :]
        print('Newmark-{} 完成'.format(chr(946)))
        x3, y3 = self.Wilson_theta()
        y3 = np.array(y3)[np.arange(0, len(y3), 10), :]
        print('Wilson-{} 完成'.format(chr(952)))
        
        with open('./results/' + self.name + '.txt', 'w') as f:
            f.write('data1 = [\n')
            for i in range(np.shape(y1)[0]):
                f.write('  [')
                for j in range(np.shape(y1)[1]):
                    f.write('[' + str(x1[j]) + ', ' + str(y1[i][j])+'], ')
                f.write('],\n')
            f.write('];\n\n')
            
            f.write('data2 = [\n')
            for i in range(np.shape(y2)[0]):
                f.write('  [')
                for j in range(np.shape(y2)[1]):
                    f.write('[' + str(x2[j]) + ', ' + str(y2[i][j])+'], ')
                f.write('],\n')
            f.write('];\n\n')

            f.write('data3 = [\n')
            for i in range(np.shape(y3)[0]):
                f.write('  [')
                for j in range(np.shape(y3)[1]):
                    f.write('[' + str(x3[j]) + ', ' + str(y3[i][j])+'], ')
                f.write('],\n')
            f.write('];\n\n')
            
        
        dataFile = ''
        with open('./results/' + self.name + '.txt', 'r') as f:
            dataFile = f.read()
            
        jsFile = ''
        with open('./js/src.js', 'r') as f:
            jsFile = f.read()
        
        jsIndexFile = dataFile + jsFile
        with open('./js/index.js', 'w') as f:
            f.write(jsIndexFile)
        
        HTMLname = './index.html'
        webbrowser.open_new_tab(HTMLname)