from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import sys
import tkinter
from tkinter import ttk
import ctypes

import sv_ttk
import darkdetect
import pywinstyles
from ujson import load, dump

from complex_mapping import *


CIRCLES = 'Окружности (Circles)'
LINES = 'Линии (Lines)'
PARABOLAS = 'Параболы (Parabolas)'

ORDER_N = 'order_n'
ZETA_0 = 'zeta_0'
KAPPA_0 = 'kappa_0'

MAX_RAD = 'max_rad'
CIRCLE_NUM = 'circle_num'
SHIFT = 'shift'

TOL = 'tol'
DISCREPANCY = 'discrepancy'


def apply_theme_to_titlebar(root):
    version = sys.getwindowsversion()

    if version.major == 10 and version.build >= 22000:
        pywinstyles.change_header_color(root, '#1c1c1c' if sv_ttk.get_theme() == 'dark' else '#fafafa')
    elif version.major == 10:
        pywinstyles.apply_style(root, 'dark' if sv_ttk.get_theme() == 'dark' else 'normal')
        root.wm_attributes('-alpha', 0.99)
        root.wm_attributes('-alpha', 1)

def config_save():
    data = {
        'model': {
            ORDER_N: int(modelframe.nametowidget(ORDER_N).get()),
            ZETA_0: modelframe.nametowidget(ZETA_0).get(),
            KAPPA_0: modelframe.nametowidget(KAPPA_0).get()
            },
        'param': {
            MAX_RAD: float(paramframe.nametowidget(MAX_RAD).get()),
            CIRCLE_NUM: int(paramframe.nametowidget(CIRCLE_NUM).get()),
            SHIFT: int(paramframe.nametowidget(SHIFT).get())
            },
        'algorithm': {
            TOL: float(algorithmframe.nametowidget(TOL).get()),
            DISCREPANCY: root.globalgetvar(DISCREPANCY)
            }
        }
    with open('config.json', 'w') as f:
        dump(data, f, indent=4)

def config_load():
    with open('config.json') as f:
        data = load(f)
        
        ordern = modelframe.nametowidget(ORDER_N)
        ordern.delete(0, 'end')
        ordern.insert(0, data['model'][ORDER_N])
        
        zeta0 = modelframe.nametowidget(ZETA_0)
        zeta0.delete(0, 'end')
        zeta0.insert(0, data['model'][ZETA_0])

        kappa0 = modelframe.nametowidget(KAPPA_0)
        kappa0.delete(0, 'end')
        kappa0.insert(0, data['model'][KAPPA_0])

        maxrad = paramframe.nametowidget(MAX_RAD)
        maxrad.delete(0, 'end')
        maxrad.insert(0, data['param'][MAX_RAD])

        circlenum = paramframe.nametowidget(CIRCLE_NUM)
        circlenum.delete(0, 'end')
        circlenum.insert(0, data['param'][CIRCLE_NUM])

        shift = paramframe.nametowidget(SHIFT)
        shift.delete(0, 'end')
        shift.insert(0, data['param'][SHIFT])

        tol = algorithmframe.nametowidget(TOL)
        tol.delete(0, 'end')
        tol.insert(0, data['algorithm'][TOL])

        discrepancy.set(data['algorithm'][DISCREPANCY])
        
def geometry_params(active_state):
    frame = paramframe
    for widget in frame.winfo_children():
        widget.destroy()
    states = {
        CIRCLES: (
            'Максимальный радиус:', MAX_RAD,
            'Количество колец:', CIRCLE_NUM,
            'Сдвиг фазы (в градусах):', SHIFT
        ),
        LINES: (
            'TEST1', 'test1'
        ),
        PARABOLAS: (
            'TEST2', 'test2'
        ),
    }
    
    state = states[active_state]
    frame.columnconfigure(0, weight=1)

    for i in range(len(state) // 2):
        
        frame.rowconfigure(i*2, weight=0)
        label1 = ttk.Label(frame, text=state[i*2])
        label1.grid(column=0, row=i*2, sticky='nsew')
        
        frame.rowconfigure(i*2+1, weight=0)
        entry1 = ttk.Entry(frame, name=state[i*2+1])
        entry1.grid(column=0, row=i*2+1, sticky='nsew')


def selected(event):
    state = event.widget.get()
    geometry_params(state)

def add_to_note(notebook, name, fig):
    frame = ttk.Frame(notebook)
    notebook.add(frame, text=name)

    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tkinter.BOTH, expand=True)

    toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=False)
    toolbar.update()
    toolbar.grid(row=2, column=1, sticky='nsew')

def start(n, zeta0, kappa0, combobox_geom, frame, notebook, progress_value):
    print('start')
    # progress_value.set(progress_value.get() + 1)
    zeta0=complex(zeta0)
    kappa0=complex(kappa0)
    cm = ComplexMapping(n, 500)
    cm.zeta0 = mp.mpc(zeta0.real, zeta0.imag)
    cm.kappa0 = mp.mpc(kappa0.real, kappa0.imag)
    if combobox_geom == CIRCLES:
        max_rad = float(frame.nametowidget(MAX_RAD).get())
        circle_num = int(frame.nametowidget(CIRCLE_NUM).get())
        shift = int(frame.nametowidget(SHIFT).get())
        cm.circles(max_rad, circle_num, shift)
        add_to_note(notebook, combobox_geom,
                    cm.plot_gradient(cm.zetaMatrix, draw_arrows=False, labels={'X': 'Real', 'Y': 'Imag'}))
    #cm.calcKappa()
    #add_to_note(notebook, cm.plot_gradient(cm.zetaMatrix, draw_arrows=False, labels={'X': 'Real', 'Y': 'Imag'}))



root = tkinter.Tk()
if hasattr(ctypes, 'windll'):
    ctypes#.windll.shcore.SetProcessDpiAwareness(1)
root.title('Complex mapping')
root.resizable(True, True)
root.geometry('1000x750')
root.eval('tk::PlaceWindow . center')

root.columnconfigure(0, weight=1)
root.columnconfigure(1, weight=4)
root.rowconfigure(0, weight=0)
root.rowconfigure(1, weight=1)
root.rowconfigure(2, weight=0)

modelframe = ttk.Frame(root, borderwidth=1, padding=10)#, border=1, relief='groove')
geometryframe = ttk.Frame(modelframe, borderwidth=1)#, padding=10, border=1, relief='solid')
algorithmframe = ttk.Frame(root, borderwidth=1, padding=10)#, border=1, relief='solid')

modelframe.grid(row=0, column=0, rowspan=3, sticky='nsew')
algorithmframe.grid(row=0, column=1, sticky='nsew')

notebook = ttk.Notebook(root)
notebook.grid(row=1, column=1, sticky='nsew')

# model frame
modelframe.columnconfigure(0, weight=1)
modelframe.columnconfigure(1, weight=1)
modelframe.rowconfigure(0, weight=0)
modelframe.rowconfigure(1, weight=0)
modelframe.rowconfigure(2, weight=0)
modelframe.rowconfigure(3, weight=0)
modelframe.rowconfigure(4, weight=0)
modelframe.rowconfigure(5, weight=0)
modelframe.rowconfigure(6, weight=0)
modelframe.rowconfigure(7, weight=1)
modelframe.rowconfigure(8, weight=0)
modelframe.rowconfigure(9, weight=0)

label0 = ttk.Label(modelframe, text=' МОДЕЛЬ', background='purple')
label0.grid(column=0, row=0, columnspan=2, sticky = 'nsew')

label1 = ttk.Label(modelframe, text = 'Порядок n:')
label1.grid(column=0, row=1, columnspan=2, sticky='nsew')
ent_n = ttk.Entry(modelframe, name=ORDER_N)
ent_n.grid(column=0, row=2, columnspan=2, sticky='nsew', pady=2)

label2 = ttk.Label(modelframe, text='Начальная zeta_0:')
label2.grid(column=0, row=3, columnspan=2, sticky='nsew')
ent_zeta_0 = ttk.Entry(modelframe, name=ZETA_0)
ent_zeta_0.grid(column=0, row=4, columnspan=2, sticky='nsew', pady=2)

label3 = ttk.Label(modelframe, text='Начальная kappa_0:')
label3.grid(column=0, row=5, columnspan=2, sticky='nsew')
ent_kappa_0 = ttk.Entry(modelframe, name=KAPPA_0)
ent_kappa_0.grid(column=0, row=6, columnspan=2, sticky='nsew', pady=2)

# geometry frame (inside model frame)
geometryframe.grid(column=0, row=7, columnspan=2, sticky='nsew', pady=20)

geometryframe.columnconfigure(0, weight=1)
geometryframe.rowconfigure(0, weight=0)
geometryframe.rowconfigure(1, weight=0)
geometryframe.rowconfigure(2, weight=1)

label4 = ttk.Label(geometryframe, text=' ГЕОМЕТРИЯ СЕТКИ', background='purple')
label4.grid(column=0, row=0, sticky='nsew')

geometries = [CIRCLES, LINES, PARABOLAS]
combobox_geom = ttk.Combobox(geometryframe, values=geometries, state='readonly')
combobox_geom.grid(column=0, row=1, sticky='nsew', pady=2)
combobox_geom.bind('<<ComboboxSelected>>', lambda event: selected(event))

#bt_start = ttk.Button(modelframe, text='Start', style='Accent.TButton', command=start(ent_n.get(), ent_zeta_0.get(), ent_kappa_0.get(), combobox_geom.get(), params, notebook))
#bt_start.grid(column=0, row=8, sticky='nsew', pady=2)

paramframe = ttk.Frame(geometryframe)#, padding=10)
paramframe.grid(column=0, row=2, sticky='nsew')

# algorithm frame
algorithmframe.columnconfigure(0, weight=0)
algorithmframe.columnconfigure(1, weight=1)
algorithmframe.columnconfigure(2, weight=0)
algorithmframe.columnconfigure(3, weight=1)
algorithmframe.columnconfigure(4, weight=0)
algorithmframe.columnconfigure(5, weight=0)
algorithmframe.columnconfigure(6, weight=0)
algorithmframe.rowconfigure(0, weight=0)
algorithmframe.rowconfigure(1, weight=0)
algorithmframe.rowconfigure(2, weight=0)

label5 = ttk.Label(algorithmframe, text=' ПАРАМЕТРЫ АЛГОРИТМА И ТОЧНОСТИ', background='purple')
label5.grid(column=0, row=0, columnspan=7, sticky='nsew')

label6 = ttk.Label(algorithmframe, text='Метод ОДУ:')
label6.grid(column=0, row=1, rowspan=2, sticky='nsew')
ody = ['Runge-Kutta (RK4)']
combobox_ody = ttk.Combobox(algorithmframe, values=ody, state='readonly')
combobox_ody.grid(column=1, row=1, rowspan=2, sticky='nsew', pady=2, padx=[5, 15])

label7 = ttk.Label(algorithmframe, text='Допуск (Tol):')
label7.grid(column=2, row=1, rowspan=2, sticky='nsew')
ent_tol = ttk.Entry(algorithmframe, name=TOL)
ent_tol.grid(column=3, row=1, rowspan=2, sticky='nsew', pady=2, padx=[5, 15])

label8 = ttk.Label(algorithmframe, text='Невязка:')
label8.grid(column=4, row=1, rowspan=2, sticky='nsew')
discrepancy = tkinter.StringVar(value='Abs', name=DISCREPANCY)
rb1 = ttk.Radiobutton(algorithmframe, text='Abs', variable=discrepancy, value='Abs')
rb1.grid(column=5, row=2, sticky='nsew', pady=2)
rb2 = ttk.Radiobutton(algorithmframe, text='Rel', variable=discrepancy, value='Rel')
rb2.grid(column=6, row=2, sticky='nsew', pady=2)


bt_save = ttk.Button(modelframe, text='Save config', style='Accent.TButton', command=config_save)
bt_save.grid(column=0, row=8, sticky='nsew', pady=[0, 5], padx=[0, 5])

bt_load = ttk.Button(modelframe, text='Load config', style='Accent.TButton', command=config_load)
bt_load.grid(column=1, row=8, sticky='nsew', pady=[0, 5], padx=[5, 0])

bt_start = ttk.Button(modelframe, text='Start', style='Accent.TButton',
                      command=lambda: start(ent_n.get(), ent_zeta_0.get(), ent_kappa_0.get(),
                                            combobox_geom.get(), paramframe, notebook, progress_value))
bt_start.grid(column=0, row=9, columnspan=2, sticky='nsew', pady=[5, 2])

progress_value = tkinter.IntVar(value=0)
#progress = ttk.Progressbar(modelframe, orient="horizontal", variable=progress_value)
#progress.grid(column=0, row=9, sticky='nsew', pady=2)


sv_ttk.set_theme('dark')
if sys.platform.startswith('win'):
    apply_theme_to_titlebar(root)
root.mainloop()
