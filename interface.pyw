from complex_mapping import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import sys
import tkinter
from tkinter import ttk

import sv_ttk
import darkdetect
import pywinstyles

def apply_theme_to_titlebar(root):
    version = sys.getwindowsversion()

    if version.major == 10 and version.build >= 22000:
        pywinstyles.change_header_color(root, "#1c1c1c" if sv_ttk.get_theme() == "dark" else "#fafafa")
    elif version.major == 10:
        pywinstyles.apply_style(root, "dark" if sv_ttk.get_theme() == "dark" else "normal")
        root.wm_attributes("-alpha", 0.99)
        root.wm_attributes("-alpha", 1)


def circles(frame):
    for widget in frame.winfo_children():
        widget.destroy()
    frame.columnconfigure(0, weight=100)
    frame.rowconfigure(0, weight=0)
    frame.rowconfigure(1, weight=0)
    frame.rowconfigure(2, weight=0)
    frame.rowconfigure(3, weight=0)
    frame.rowconfigure(4, weight=0)
    frame.rowconfigure(5, weight=0)

    label1 = ttk.Label(frame, text="Макс. радиус:")
    label1.grid(column=0, row=0, sticky='nsew')
    label2 = ttk.Label(frame, text="Количество колец:")
    label2.grid(column=0, row=2, sticky='nsew')
    label3 = ttk.Label(frame, text="Сдвиг фазы (град):")
    label3.grid(column=0, row=4, sticky='nsew')

    max_rad = ttk.Entry(frame, name="max_rad")
    max_rad.grid(column=0, row=1, sticky='nsew')
    circl_num = ttk.Entry(frame, name="circl_num")
    circl_num.grid(column=0, row=3, sticky='nsew')
    shift = ttk.Entry(frame, name="shift")
    shift.grid(column=0, row=5, sticky='nsew')

def selected(event, frame):
    state = event.widget.get()
    if (state == "Окружности (Circles)"):
        circles(frame)

def add_to_note(notebook, fig):
    frame = ttk.Frame(notebook)
    notebook.add(frame)
    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tkinter.BOTH, expand=True)

def start(n, zeta0, kappa0, combobox_geom, frame, notebook):
    print("start")
    zeta0=complex(zeta0)
    kappa0=complex(kappa0)
    cm = ComplexMapping(n, 500)
    cm.zeta0 = mp.mpc(zeta0.real, zeta0.imag)
    cm.kappa0 = mp.mpc(kappa0.real, kappa0.imag)
    if combobox_geom == "Окружности (Circles)":
        max_rad = float(frame.nametowidget("max_rad").get())
        circl_num = int(frame.nametowidget("circl_num").get())
        shift = int(frame.nametowidget("shift").get())
        cm.circles(max_rad, circl_num, shift)
        add_to_note(notebook, cm.plot_gradient(cm.zetaMatrix, draw_arrows=False, labels={'X': 'Real', 'Y': 'Imag'}))
    #cm.calcKappa()
    #add_to_note(notebook, cm.plot_gradient(cm.zetaMatrix, draw_arrows=False, labels={'X': 'Real', 'Y': 'Imag'}))



def main():
    root = tkinter.Tk()
    root.title('Test app')
    root.resizable(True, True)
    root.geometry('1000x600')
    root.eval('tk::PlaceWindow . center')

    root.columnconfigure(0, weight=20)
    root.columnconfigure(1, weight=80)
    root.rowconfigure(0, weight=0)
    root.rowconfigure(1, weight=100)

    modelframe = ttk.Frame(root, borderwidth = 1, padding = 10, border = 1, relief = 'solid')
    geometryframe = ttk.Frame(modelframe, borderwidth=1, padding=10, border = 1, relief = 'solid')
    algorithmframe = ttk.Frame(root, borderwidth=1, padding=10, border = 1, relief='solid')
    modelframe.grid(row=0, column=0, rowspan=2, sticky = 'nsew')
    algorithmframe.grid(row=0, column=1, sticky = 'nsew')

    notebook = ttk.Notebook(root)
    notebook.grid(row=1, column=1, sticky='nsew')

    #modelframe
    modelframe.columnconfigure(0, weight=100)
    modelframe.rowconfigure(0, weight=0)
    modelframe.rowconfigure(1, weight=0)
    modelframe.rowconfigure(2, weight=0)
    modelframe.rowconfigure(3, weight=0)
    modelframe.rowconfigure(4, weight=0)
    modelframe.rowconfigure(5, weight=0)
    modelframe.rowconfigure(6, weight=0)
    modelframe.rowconfigure(7, weight=80)
    modelframe.rowconfigure(8, weight=20)

    label0 = ttk.Label(modelframe, text="МОДЕЛЬ", background="purple")
    label0.grid(column=0, row=0, sticky = 'nsew')
    label1 = ttk.Label(modelframe, text = "Порядок n:")
    label1.grid(column=0, row=1, sticky='nsew')
    label2 = ttk.Label(modelframe, text="Начальная zeta_0:")
    label2.grid(column=0, row=3, sticky='nsew')
    label3 = ttk.Label(modelframe, text="Начальная kappa_0:")
    label3.grid(column=0, row=5, sticky='nsew')

    ent_n = ttk.Entry(modelframe)
    ent_n.grid(column=0, row=2, sticky='nsew', pady=2)
    ent_zeta_0 = ttk.Entry(modelframe)
    ent_zeta_0.grid(column=0, row=4, sticky='nsew', pady=2)
    ent_kappa_0 = ttk.Entry(modelframe)
    ent_kappa_0.grid(column=0, row=6, sticky='nsew', pady=2)

    geometryframe.grid(column=0, row=7, sticky='nsew', pady=2)



    geometryframe.columnconfigure(0, weight=100)
    geometryframe.rowconfigure(0, weight=1)
    geometryframe.rowconfigure(1, weight=1)
    geometryframe.rowconfigure(2, weight=98)

    label4 = ttk.Label(geometryframe, text="ГЕОМЕТРИЯ СЕТКИ", background="purple")
    label4.grid(column=0, row=0, sticky='nsew')

    geometries = ["Окружности (Circles)", "Линии (Lines)", "Параболы (Parabols)"]
    combobox_geom = ttk.Combobox(geometryframe, values=geometries, state="readonly")
    combobox_geom.grid(column=0, row=1, sticky='nsew', pady=2)
    combobox_geom.bind("<<ComboboxSelected>>", lambda event: selected(event, params_gem))

    #bt_start = ttk.Button(modelframe, text='Start', style='Accent.TButton', command=start(ent_n.get(), ent_zeta_0.get(), ent_kappa_0.get(), combobox_geom.get(), params, notebook))
    #bt_start.grid(column=0, row=8, sticky='nsew', pady=2)

    params_gem = ttk.Frame(geometryframe, padding=10)
    params_gem.grid(column=0, row=2, sticky='nsew')

    #algorithmframe
    algorithmframe.rowconfigure(0, weight=0)
    algorithmframe.rowconfigure(1, weight=0)
    algorithmframe.rowconfigure(2, weight=0)
    algorithmframe.columnconfigure(0, weight=0)
    algorithmframe.columnconfigure(1, weight=45)
    algorithmframe.columnconfigure(2, weight=0)
    algorithmframe.columnconfigure(3, weight=45)
    algorithmframe.columnconfigure(4, weight=0)
    algorithmframe.columnconfigure(5, weight=10)

    label5 = ttk.Label(algorithmframe, text="ПАРАМЕТРЫ АЛГОРИТМА И ТОЧНОСТИ", background="purple")
    label5.grid(column=0, row=0, columnspan=6, sticky='nsew')
    label6 = ttk.Label(algorithmframe, text="Метод ОДУ:")
    label6.grid(column=0, row=1, rowspan=2, sticky='nsew')
    label7 = ttk.Label(algorithmframe, text="Допуск (Tol):")
    label7.grid(column=2, row=1, rowspan=2, sticky='nsew')
    label8 = ttk.Label(algorithmframe, text="Невязка:")
    label8.grid(column=4, row=1, rowspan=2, sticky='nsew')

    ody = ["RK4", "???"]
    combobox_ody = ttk.Combobox(algorithmframe, values=ody, state="readonly")
    combobox_ody.grid(column=1, row=1, rowspan=2, sticky='nsew', pady=2)

    ent_tol = ttk.Entry(algorithmframe)
    ent_tol.grid(column=3, row=1, rowspan=2, sticky='nsew', pady=2)

    nevyazka = tkinter.StringVar(value="Abs")
    rb1 = ttk.Radiobutton(algorithmframe, text="Abs", variable=nevyazka, value="Abs")
    rb1.grid(column=5, row=2, sticky='nsew', pady=2)
    rb2 = ttk.Radiobutton(algorithmframe, text="Rel", variable=nevyazka, value="Rel")
    rb2.grid(column=5, row=3, sticky='nsew', pady=2)

    bt_start = ttk.Button(modelframe, text='Start', style='Accent.TButton', command=lambda: start(ent_n.get(), ent_zeta_0.get(), ent_kappa_0.get(), combobox_geom.get(), params_gem, notebook))
    bt_start.grid(column=0, row=8, sticky='nsew', pady=2)




    '''
    sv_ttk.set_theme('light')
    if sys.platform.startswith('win'):
        apply_theme_to_titlebar(root)
    '''
    root.mainloop()


if __name__ == '__main__':
    main()

