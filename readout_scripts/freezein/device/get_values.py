import json
import numpy as np


def decode_value(val):
    dt_coarse = 1/175*1e3 # in ns
    dt_fine = 1/350/8*1e3
    if val[:2]=='0x':
        value = bin(int(val[2:], base=16))[2:].zfill(40)
    channel = value[-37-1:-35-1]
    edge_bit = value[-33-1]
    #
    fine_time = value[-2-1:]
    phase_bit = value[-3-1]
    coarse_time = value[-32-1:-4]
    epoch_rollover = value[-35-1]
    #
    decode_error = value[-34-1]

    time_ns = int(coarse_time, base=2)*dt_coarse + int(fine_time, base=2)*dt_fine + dt_coarse/2 * int(phase_bit=='1')
    return (int(channel,base=2), int(edge_bit), time_ns, int(epoch_rollover))


def get_timestamps(in_):
    dtype=[('ch', int),
       ('edge', int),
       ('time_ns', float),
       ('epoch_ro', int),
       ('epoch', int)
      ]

    out_ = np.zeros(len(in_), dtype=dtype)
    for i_ in range(len(in_)):
        out_[i_] = decode_value(in_[i_])+(0,)
    out_['epoch'] = np.cumsum(out_['epoch_ro'])
    out_['time_ns'] += (2**29/0.175 ) * out_['epoch']
    return out_[out_['epoch_ro']==0]


def getADCreadings(in_):
    dtype=[('base', int),
       ('peak', int),
       ('sum', int), 
      ]
    out_ = np.zeros(len(in_), dtype=dtype)
    for i_ in range(len(in_)):
        out_['base'][i_] = int( '0x' + (in_[i_].encode("utf-8")).decode()[-4:], base=16)
        out_['peak'][i_] = int( '0x' + (in_[i_].encode("utf-8")).decode()[-8:-4], base=16)
        out_['sum'][i_]  = int( '0x' + (in_[i_].encode("utf-8")).decode()[-16:-8], base=16)
    return out_


def convert_to_perTrig(indata, ntrigs = None, timerange = (-200, 100000)):
    ntrigs = int(indata['tdc3'].shape[0]/2)
    outdata = {}
    outdata['tdc0'] = []
    outdata['tdc2'] = []
    outdata['tdc3'] = np.zeros( (ntrigs, 2), dtype = indata['tdc3'].dtype )
    #print(outdata['tdc3'])
    for i_ in range(0, ntrigs):
        outdata['tdc3'][i_] = indata['tdc3'][2*i_:2*i_+2]
        #print(outdata['tdc3'][i_])
        for fifo_ in ['tdc0', 'tdc2']:
            selbool_ = ( ( (indata[fifo_]['time_ns'] - outdata['tdc3']['time_ns'][i_, 0] ) > timerange[0] )*
                         ( (indata[fifo_]['time_ns'] - outdata['tdc3']['time_ns'][i_, 1] ) < timerange[1] )
                       )
            outdata[fifo_].append(np.array(indata[fifo_][selbool_]))
            del selbool_
    for fifo_ in ['adcA','adcB']:
        outdata[fifo_]= np.array(indata[fifo_])
    return(outdata)


def get_SiPM_dt(data):
    out_ = {}
    for fifo_  in ['tdc0', 'tdc2', 'tdc3']:
        out_[fifo_+"_dt"] = np.zeros(len(data[fifo_]))
        for i_ in range(len(out_[fifo_+"_dt"])):
            rise_bool = data[fifo_][i_]['edge'] == 1
            fall_bool = data[fifo_][i_]['edge'] == 0
            out_[fifo_+"_dt"][i_] = (data[fifo_][i_]['time_ns'][fall_bool][0] - data[fifo_][i_]['time_ns'][rise_bool][0])
    return(out_)



def extract(dataset, side):
    data = {}
    if side == 'm':
        val = 'vals_m'
    else:
        val = 'vals_s'
    for fifo_ in ['tdc0','tdc2','tdc3']:
        data[fifo_] = get_timestamps(dataset[fifo_][val])
    for fifo_ in ['adcA', 'adcB']:
        data[fifo_] = getADCreadings(dataset[fifo_][val])
    data = convert_to_perTrig(data)
    data.update(get_SiPM_dt(data))
    return(data)


def extract_dark(dataset, side):
    data = {}
    if side == 'm':
        val = 'vals_m'
    else:
        val = 'vals_s'
    for fifo_ in ['adcA', 'adcB']:
        data[fifo_] = np.array(getADCreadings(dataset[fifo_][val]))
    return(data)
