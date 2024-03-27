from appJar import gui
import os

BASEPATH = "src/scripts/methylation"


def pick_script(path):

    def press(name):
        next_path = f"{path}/{name}"
        app.stop()
        if os.path.isdir(next_path):
            pick_script(next_path)
        else:
            with open(next_path) as script:
                print(next_path)
                exec(script.read())

    files = [filename for filename in os.listdir(path) if filename[0] != "_"]
    app = gui()
    app.addOptionBox("Script", files)
    def run():
        filename = app.getOptionBox("Script")
        press(filename)
    app.addButton("run", run)
    app.go()

pick_script(BASEPATH)

