
import subprocess as sb
import win32com.client
import pythoncom


def make_SW_part(Spring):
    startSW()
    sw = connectToSW()
    newFile(sw, "newPart.SLDPRT")


def startSW():
    ## Starts Solidworks
    SW_PROCESS_NAME = r'C:/Program Files/SOLIDWORKS Corp/SOLIDWORKS/SLDWORKS.exe'
    sb.Popen(SW_PROCESS_NAME)

def shutSW():
    ## Kills Solidworks
    sb.call('Taskkill /IM SLDWORKS.exe /F')

def connectToSW():
    ## With Solidworks window open, connects to application
    sw = win32com.client.Dispatch("SLDWORKS.Application")
    return sw

def openFile(sw, Path):
    ## With connection established (sw), opens part, assembly, or drawing file
    f = sw.getopendocspec(Path)
    model = sw.opendoc7(f)
    return model

def newFile(sw, Path):
    template = "C:\ProgramData\SolidWorks\SOLIDWORKS 2014\templates\Part.prtdot"
    model = sw.NewDocument(template, 0,0,0)
    model.SaveAs(Path, 0, 2)
    return model

def updatePrt(model):

    ## Rebuilds the active part, assembly, or drawing (model)
    model.EditRebuild3
