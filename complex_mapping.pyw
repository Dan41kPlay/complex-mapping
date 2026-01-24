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


def main():
    root = tkinter.Tk()
    root.title('Test app')
    root.resizable(False, False)
    root.geometry('500x300')
    root.eval('tk::PlaceWindow . center')

    b1 = ttk.Button(root, text='Test 1', style='Accent.TButton')
    b1.place(x=175, y=75, width=150, height=50)

    b2 = ttk.Button(root, text='Test 2')
    b2.place(x=175, y=175, width=150, height=50)

    sv_ttk.set_theme('dark')
    if sys.platform.startswith('win'):
        apply_theme_to_titlebar(root)

    root.mainloop()


if __name__ == '__main__':
    main()
