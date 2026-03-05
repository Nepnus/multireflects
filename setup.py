import os
import sys
import shutil
import sysconfig
import subprocess

try:
    import numpy
except:
    raise Exception("Please install numpy.")
try:
    import scipy
except:
    raise Exception("Please install scipy.")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
pythonlib_path = sysconfig.get_path('purelib')
install_path = os.path.join(pythonlib_path, "multireflects")

if sys.platform.startswith("linux"):
    try:
        subprocess.run("gcc --version", check=True, shell=True, stdout=subprocess.DEVNULL)
    except:
        raise Exception("Please install and configure the gcc compiler correctly.")
    try:
        path1 = os.path.join(BASE_DIR, "CubicSplineInterp.c")
        path2 = os.path.join(BASE_DIR, "libCubicSplineInterp.so")
        subprocess.run(f"gcc -O2 -shared -fPIC -Wl,-rpath,\'$ORIGIN\' -I {BASE_DIR} {path1} -o {path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "Fin_interp.c")
        path2 = os.path.join(BASE_DIR, "libFin_interp.so")
        subprocess.run(f"gcc -O2 -shared -fPIC -Wl,-rpath,\'$ORIGIN\' -I {BASE_DIR} -L {BASE_DIR} -lCubicSplineInterp {path1} -o {path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "Units_core.c")
        path2 = os.path.join(BASE_DIR, "Units_core.so")
        subprocess.run(f"gcc -O2 -shared -fPIC -Wl,-rpath,\'$ORIGIN\' {path1} -o {path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "Fin_core.c")
        path2 = os.path.join(BASE_DIR, "Fin_core.so")
        subprocess.run(f"gcc -O2 -shared -fPIC -Wl,-rpath,\'$ORIGIN\' {path1} -o {path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "lc_mode0.c")
        path2 = os.path.join(BASE_DIR, "lc_mode0.so")
        subprocess.run(f"gcc -O2 -shared -fPIC -Wl,-rpath,\'$ORIGIN\' -I {BASE_DIR} -L {BASE_DIR} -lCubicSplineInterp -lFin_interp {path1} -o {path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "lc_mode1.c")
        path2 = os.path.join(BASE_DIR, "lc_mode1.so")
        subprocess.run(f"gcc -O2 -shared -fPIC -Wl,-rpath,\'$ORIGIN\' -I {BASE_DIR} -L {BASE_DIR} -lCubicSplineInterp {path1} -o {path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
    except:
        raise Exception("Compilation failed.")
    if not os.path.exists(install_path):
        os.makedirs(install_path, exist_ok=True)
    filenamelist = os.listdir(BASE_DIR)
    for filename in filenamelist:
        if filename == "setup.py" or filename == "check.py":
            continue
        if ".so" not in filename and ".py" not in filename:
            continue
        path1 = os.path.join(BASE_DIR, filename)
        path2 = os.path.join(install_path, filename)
        shutil.copy2(path1, path2)
    
elif sys.platform.startswith("darwin"):
    try:
        subprocess.run("clang --version", check=True, shell=True, stdout=subprocess.DEVNULL)
    except:
        raise Exception("Please install and configure the clang compiler correctly.")
    try:
        path1 = os.path.join(BASE_DIR, "CubicSplineInterp.c")
        path2 = os.path.join(BASE_DIR, "libCubicSplineInterp.dylib")
        subprocess.run(f"clang -O2 -dynamiclib -fPIC -arch arm64 -arch x86_64 -Wl,-rpath,@loader_path -install_name @rpath/libCubicSplineInterp.dylib -I {BASE_DIR} {path1} -o {path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "Fin_interp.c")
        path2 = os.path.join(BASE_DIR, "libFin_interp.dylib")
        subprocess.run(f"clang -O2 -dynamiclib -fPIC -arch arm64 -arch x86_64 -Wl,-rpath,@loader_path -install_name @rpath/libFin_interp.dylib -I {BASE_DIR} -L {BASE_DIR} -lCubicSplineInterp {path1} -o {path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "Units_core.c")
        path2 = os.path.join(BASE_DIR, "Units_core.dylib")
        subprocess.run(f"clang -O2 -dynamiclib -fPIC -arch arm64 -arch x86_64 -Wl,-rpath,@loader_path -install_name @rpath/Units_core.dylib {path1} -o {path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "Fin_core.c")
        path2 = os.path.join(BASE_DIR, "Fin_core.dylib")
        subprocess.run(f"clang -O2 -dynamiclib -fPIC -arch arm64 -arch x86_64 -Wl,-rpath,@loader_path -install_name @rpath/Fin_core.dylib {path1} -o {path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "lc_mode0.c")
        path2 = os.path.join(BASE_DIR, "lc_mode0.dylib")
        subprocess.run(f"clang -O2 -dynamiclib -fPIC -arch arm64 -arch x86_64 -Wl,-rpath,@loader_path -install_name @rpath/lc_mode0.dylib -I {BASE_DIR} -L {BASE_DIR} -lCubicSplineInterp -lFin_interp {path1} -o {path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "lc_mode1.c")
        path2 = os.path.join(BASE_DIR, "lc_mode1.dylib")
        subprocess.run(f"clang -O2 -dynamiclib -fPIC -arch arm64 -arch x86_64 -Wl,-rpath,@loader_path -install_name @rpath/lc_mode1.dylib -I {BASE_DIR} -L {BASE_DIR} -lCubicSplineInterp {path1} -o {path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
    except:
        raise Exception("Compilation failed.")
    if not os.path.exists(install_path):
        os.makedirs(install_path, exist_ok=True)
    filenamelist = os.listdir(BASE_DIR)
    for filename in filenamelist:
        if filename == "setup.py" or filename == "check.py":
            continue
        if ".dylib" not in filename and ".py" not in filename:
            continue
        path1 = os.path.join(BASE_DIR, filename)
        path2 = os.path.join(install_path, filename)
        shutil.copy2(path1, path2)

elif sys.platform.startswith("win"):
    try:
        subprocess.run("cl", check=True, shell=True, stdout=subprocess.DEVNULL)
    except:
        raise Exception("Please install and configure the msvc compiler correctly.")
    try:
        path1 = os.path.join(BASE_DIR, "CubicSplineInterp.c")
        path2 = os.path.join(BASE_DIR, "libCubicSplineInterp.dll")
        subprocess.run(f"cl /O2 /LD /I {BASE_DIR} {path1} /link /OUT:{path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "Fin_interp.c")
        path2 = os.path.join(BASE_DIR, "libFin_interp.dll")
        subprocess.run(f"cl /O2 /LD /I {BASE_DIR} {path1} CubicSplineInterp.obj /link /LIBPATH:{BASE_DIR} /OUT:{path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "Units_core.c")
        path2 = os.path.join(BASE_DIR, "Units_core.dll")
        subprocess.run(f"cl /O2 /LD {path1} /link /OUT:{path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "Fin_core.c")
        path2 = os.path.join(BASE_DIR, "Fin_core.dll")
        subprocess.run(f"cl /O2 /LD {path1} /link /OUT:{path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "lc_mode0.c")
        path2 = os.path.join(BASE_DIR, "lc_mode0.dll")
        subprocess.run(f"cl /O2 /LD /I {BASE_DIR} {path1} CubicSplineInterp.obj Fin_interp.obj /link /LIBPATH:{BASE_DIR} /OUT:{path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
        path1 = os.path.join(BASE_DIR, "lc_mode1.c")
        path2 = os.path.join(BASE_DIR, "lc_mode1.dll")
        subprocess.run(f"cl /O2 /LD /I {BASE_DIR} {path1} CubicSplineInterp.obj /link /LIBPATH:{BASE_DIR} /OUT:{path2}", 
                       check=True, shell=True, stdout=subprocess.DEVNULL)
    except:
        raise Exception("Compilation failed.")
    if not os.path.exists(install_path):
        os.makedirs(install_path, exist_ok=True)
    filenamelist = os.listdir(BASE_DIR)
    for filename in filenamelist:
        if filename == "setup.py" or filename == "check.py":
            continue
        if ".dll" not in filename and ".py" not in filename:
            continue
        path1 = os.path.join(BASE_DIR, filename)
        path2 = os.path.join(install_path, filename)
        shutil.copy2(path1, path2)

else:
    raise Exception("Unsupportted system.")

print("multireflects version 1.0 is installed.")



