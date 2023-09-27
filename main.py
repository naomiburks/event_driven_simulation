from appJar import gui
import os

def press(name):
    with open(f"src/scripts/{name}") as f:
        exec(f.read())


files = [fname for fname in os.listdir("src/scripts/") if fname[0] != "_"]

app = gui()

app.addLabel("title", "Pick a Script")
app.setLabelBg("title", "white")
for fname in files:
    app.addButton(fname, press)

app.go()
