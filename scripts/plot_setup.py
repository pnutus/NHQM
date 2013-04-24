import matplotlib

def plot_init(font_size =12,tick_size =8, tick_pad =6):
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : font_size}
            

    matplotlib.rc('font', **font)

    matplotlib.rc('xtick', labelsize=tick_size) 
    matplotlib.rc('ytick', labelsize=tick_size)
    
    matplotlib.rcParams['xtick.major.pad'] = tick_pad
    matplotlib.rcParams['ytick.major.pad'] = tick_pad
