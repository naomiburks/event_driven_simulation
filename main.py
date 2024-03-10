from appJar import gui
import os

def press(name):
    with open(f"src/scripts/{name}") as script:
        exec(script.read())


files = [filename for filename in os.listdir("src/scripts/") if filename[0] != "_"]

app = gui()
app.addOptionBox("Script", files)

def run():
    filename = app.getOptionBox("Script")
    press(filename)

app.addButton("run", run)
#app.addLabel("title", "Pick a Script")
#app.setLabelBg("title", "white")
#for filename in files:
#    app.addButton(filename, press)

app.go()
