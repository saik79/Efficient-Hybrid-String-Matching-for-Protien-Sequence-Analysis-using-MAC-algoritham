from tkinter import messagebox
from tkinter import *
from tkinter import simpledialog
import tkinter
from tkinter import filedialog
import matplotlib.pyplot as plt
from tkinter.filedialog import askopenfilename
import numpy as np
from CustomButton import TkinterCustomButton
import timeit
import pandas as pd
import decimal
import matplotlib.pyplot as plt
import webbrowser


main = Tk()
main.title("Efficient Hybrid String Matching Algorithm for Protein Sequence Analysis")
main.geometry("1300x1200")


global filename, dataset, user_input
global computation_time

def uploadDataset():
    global filename, dataset
    text.delete('1.0', END)
    filename = filedialog.askopenfilename(initialdir = "Dataset")
    text.insert(END,"Dataset Loaded\n\n")
    dataset= pd.read_csv(filename)
    text.insert(END,str(dataset))

def berry_ravindran(text, pattern):
    n = len(text)
    m = len(pattern)
    if m == 0:
        return list(range(n + 1))
    if m > n:
        return []
    occurrences = []
    i = 0
    while i <= n - m:
        j = 0
        while j < m:
            if i + j >= n or text[i + j] != pattern[j]:
                break  # Mismatch or end of text
            j += 1
        if j == m:  # Pattern found
            occurrences.append(i)
            i += 1  # Standard shift after a match
        else:
            # Heuristic shift based on mismatched character (complex loop)
            if j > 0 and i + j < n:
                mismatched_char = text[i + j]
                last_occurrence = -1
                for k in range(j - 1, -1, -1):
                    if pattern[k] == mismatched_char:
                        last_occurrence = k
                        break
                if last_occurrence != -1:
                    i += max(1, j - last_occurrence - 1)
                else:
                    i += max(1, j)
            else:
                i += 1 #Handles j == 0 or i+j >= n cases.
    return occurrences

def index_based_shifting(text, pattern):
    n = len(text)
    m = len(pattern)
    if m == 0:
        return list(range(n + 1))  # Empty pattern matches everywhere
    if m > n:
        return []  # Pattern longer than text, no match
    occurrences = []
    i = 0
    while i <= n - m:
        j = 0
        while j < m and text[i + j] == pattern[j]:
            j += 1
        if j == m:
            occurrences.append(i)
            i += 1 # shift by one after a match.
        else:
            i += 1 # shift by one if no match.
    return occurrences

def mac(text, pattern):
    n = len(text)
    m = len(pattern)
    if m == 0:
        return list(range(n + 1))  # Empty pattern matches everywhere
    if m > n:
        return []  # Pattern longer than text, no match
    occurrences = []
    i = 0
    while i <= n - m:
        if text[i:i + m] == pattern:  # Berry-Ravindran  comparison
            occurrences.append(i)
            i += 1  # Index-based shifting (shift by one)
        else:
            i += 1  # Index-based shifting (shift by one)
    return occurrences

def graph():
    global computation_time
    height = computation_time
    bars = ['Berry Ravindran', 'IBS', 'Hybrid MAC']
    y_pos = np.arange(len(bars))
    plt.figure(figsize = (4, 3)) 
    plt.bar(y_pos, height)
    plt.xticks(y_pos, bars)
    plt.xlabel("Algorithm Names")
    plt.ylabel("Execution Time")
    plt.title("Computation Cost Graph")
    plt.show()

def generateOutput(output_arr, index):
    output = '<table border=1 align=center>'
    output+= '<tr><th>Sequence</th><th>Pattern Matched Index</th><th>Sequence Class Name</th></tr>'
    for i in range(len(output_arr)):
        out = output_arr[i]
        sequence = out[0]
        temp = ""
        start = int(out[1])
        end = int(out[1]) + index
        for j in range(len(sequence)):
            if j >= start and j < end:
                temp += '<font size=3 color=red>'+sequence[j]+'</font>'
            else:
                temp += sequence[j]
        output+='<tr><td>'+temp+'</td><td>'+str(out[1])+'</td><td>'+str(out[2])+'</td><td></tr>'
    output+='</table></body></html>'
    f = open("output.html", "w")
    f.write(output)
    f.close()
    webbrowser.open("output.html",new=1)  

def runBerry():
    global computation_time, dataset
    text.delete('1.0', END)
    computation_time = []
    data = dataset.values
    matches = []
    pattern = tf1.get().strip()
    if len(pattern) > 0:
        start = timeit.default_timer()
        for i in range(0, 100):
            sequence = data[i,1]
            match = berry_ravindran(sequence, pattern)
            if len(match) > 0:
                matches.append([sequence, match[0], data[i,9]])
        end = timeit.default_timer()
        cost = round(decimal.Decimal(end - start), 7)
        computation_time.append(cost)
        text.insert(END,"Berry String Matching Computation Time = "+str(cost)+"\n")
        generateOutput(matches, len(pattern))
    else:
        text.insert(END,"Pattern must be entered\n")
        
def runIBS():
    global computation_time, dataset
    data = dataset.values
    matches = []
    pattern = tf1.get().strip()
    if len(pattern) > 0:
        start = timeit.default_timer()
        for i in range(0, 100):
            sequence = data[i,1]
            match = index_based_shifting(sequence, pattern)
            if len(match) > 0:
                matches.append([sequence, match[0], data[i,9]])
        end = timeit.default_timer()
        cost = round(decimal.Decimal(end - start), 7)
        computation_time.append(cost)
        text.insert(END,"IBS String Matching Computation Time = "+str(cost)+"\n")
        generateOutput(matches, len(pattern))
    else:
        text.insert(END,"Pattern must be entered\n")
    
def runMAC():
    global computation_time, dataset
    data = dataset.values
    matches = []
    pattern = tf1.get().strip()
    if len(pattern) > 0:
        start = timeit.default_timer()
        for i in range(0, 100):
            sequence = data[i,1]
            match = mac(sequence, pattern)
            if len(match) > 0:
                matches.append([sequence, match[0], data[i,9]])
        end = timeit.default_timer()
        cost = round(decimal.Decimal(end - start), 7)
        computation_time.append(cost)
        text.insert(END,"Hybrid MAC String Matching Computation Time = "+str(cost)+"\n")
        generateOutput(matches, len(pattern))
    else:
        text.insert(END,"Pattern must be entered\n")
    

font = ('times', 15, 'bold')
title = Label(main, text='Efficient Hybrid String Matching Algorithm for Protein Sequence Analysis')
title.config(bg='HotPink4', fg='yellow2')  
title.config(font=font)           
title.config(height=3, width=120)       
title.place(x=0,y=5)

font1 = ('times', 13, 'bold')
ff = ('times', 12, 'bold')

uploadButton = TkinterCustomButton(text="Load Protein Database", width=300, corner_radius=5, command=uploadDataset)
uploadButton.place(x=50,y=100)

l1 = Label(main, text='User Input:')
l1.config(font=font1)
l1.place(x=50,y=150)

tf1 = Entry(main,width=60)
tf1.config(font=font1)
tf1.place(x=200,y=150)

preprocessButton = TkinterCustomButton(text="Berry Ravindran Processing", width=270, corner_radius=5, command=runBerry)
preprocessButton.place(x=50,y=200)

fsButton = TkinterCustomButton(text="Run IBS Based Processing", width=270, corner_radius=5, command=runIBS)
fsButton.place(x=370,y=200)

splitButton = TkinterCustomButton(text="Run Hybrid MAC Processing", width=270, corner_radius=5, command=runMAC)
splitButton.place(x=690,y=200)

cnnButton = TkinterCustomButton(text="Comparison Graph", width=270, corner_radius=5, command=graph)
cnnButton.place(x=1010,y=200)

font1 = ('times', 13, 'bold')
text=Text(main,height=20,width=130)
scroll=Scrollbar(text)
text.configure(yscrollcommand=scroll.set)
text.place(x=10,y=250)
text.config(font=font1)

main.config(bg='plum2')
main.mainloop()
