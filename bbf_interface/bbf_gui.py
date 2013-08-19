from Tkinter import *
import numpy as np
master=Tk()
#var=IntVar()

Nlinac=1
Nmeas=2
Ncont=2
def cbcommand(statein,event):
    print statein

state=np.zeros((Nmeas*Nlinac,Ncont*Nlinac),dtype=type(IntVar()))
for k in range(Nmeas*Nlinac):
    for l in range(Ncont*Nlinac):
        print k
        print l
        c = Checkbutton(master, variable=state[k,l])#, command=cbcommand(state[k,l],event))
        c.grid(column=k*10, row=l*10, sticky=(W, E))
        c.pack()
    
mainloop()
