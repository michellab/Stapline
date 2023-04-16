import


def residual(vars, x, data, eps_data):
  amp = vars[0]phaseshift = vars[1]
  freq = vars[2]decay = vars[3]
  model = amp*sin(x*freq  + phaseshift)*exp(-x*x*decay)
  return(data-model)/eps_data



#from scipy.optimize importleastsqvars = [10.0, 0.2, 3.0, 0.007]out = leastsq(residual, vars, args=(x, data, eps_data))



from lmfit import minimize, Parameters

def residual(params, x, data, eps_data):
  amp = params['amp'].value
  pshift = params['phase'].value
  freq = params['frequency'].value
  decay = params['decay'].value
  model = amp*sin(x*freq  + pshift)*exp(-x*x*decay)
  return(data-model)/eps_data

params = Parameters()
params.add('amp', value=10)
params.add('decay', value=0.007)
params.add('phase', value=0.2)
params.add('frequency', value=3.0)
out = minimize(residual, params, args=(x, data, eps_data))


params = Parameters()
params.add('amp', value=10, max=30)
params.add('order', value=0.0, min=0.0, max=7.0 ,  expr='is int')
params.add('phase', value=0  , expr='in [0,pi]')
params.add('frequency', value=0.0, max=10)
params['amp'].vary = False
params['decay'].min = 0.10
