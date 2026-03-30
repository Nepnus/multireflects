## 📈 About

`multireflects` is a Python library designed for generating multi-band reflection-dominated light curves, with computational acceleration implemented in C. It can be integrated with MCMC sampling libraries, such as `emcee`, to explore the orbital parameters of binary systems dominated by reflection effects.

The main advantages of `multireflects` are:

(1) A single execution can simultaneously produce light curves in multiple passbands.

(2) For the same input parameters, its runtime performance (excluding initialization time) is at least one order of magnitude faster than commonly used reflection light-curve generators such as `phoebe`, and on average about two orders of magnitude faster.

(3) The stellar surface is divided into a fine mesh (over 10,000 surface elements), resulting in smooth light curves.

The main limitations are:

(1) The use of ellipsoidal approximations and relatively coarse surface gravity interpolation can introduce significant errors when modeling ellipsoidal-variation-dominated light curves, particularly in regions where gravity darkening is strong.

(2) Users must properly configure a C compiler and possess some basic knowledge of compilation to successfully install and use the library.

## 🚀 Installation

This project requires Python, a C compiler, and several scientific libraries.

---

### 1️⃣ Install Python Dependencies

You need to install **NumPy** and **SciPy** first:

```bash
pip install numpy
pip install scipy
```

---

### 2️⃣ Configure a C Compiler (Platform Dependent)

Since this project includes a C backend, you must have a working C compiler configured on your system.

---

## 🐧 Linux

Make sure **GCC** is installed.

Check your GCC installation:

```bash
gcc --version
```

---

## 🍎 macOS

Make sure **Clang** is installed.

Check your Clang installation:

```bash
clang --version
```

If Clang is not installed, download and install **Xcode** from the Mac App Store.  
The required compiler will be installed automatically.

---

## 🖥️ Windows

You need to configure the **MSVC** compiler.

Check whether MSVC is available:

```bash
cl
```

If the command is not recognized, follow the steps below to install and configure MSVC.

---

### 🔧 Install and Configure MSVC (Windows)

#### Step 1  
Go to:

https://visualstudio.microsoft.com/downloads/

Download **Visual Studio Installer** (⚠️ not the full Visual Studio IDE).

---

#### Step 2  
Run the installer and select:

```
Visual Studio Build Tools
  → Individual components
     → MSVC vxxx - VS xxxx C++ x86/64 build tools
     → Windows xx SDK
```

Install the selected components.

---

#### Step 3  
Locate the installation directories (usually):

```
C:\Program Files (x86)\Microsoft Visual Studio\
C:\Program Files (x86)\Windows Kits\
```

---

#### Step 4  
Add the following paths to your system **Path** environment variable:

```
(your VS path)\(xxxx)\BuildTools\VC\Tools\MSVC\(xx.xx.xxxxx)\bin\Hostx64\x64
(your WinKits path)\(xx)\bin\(xx.x.xxxxx.x)\x64
```

---

#### Step 5  
Add the following paths to the **LIB** environment variable:

```
(your VS path)\(xxxx)\BuildTools\VC\Tools\MSVC\(xx.xx.xxxxx)\lib\x64
(your WinKits path)\(xx)\Lib\(xx.x.xxxxx.x)\ucrt\x64
(your WinKits path)\(xx)\Lib\(xx.x.xxxxx.x)\um\x64
```

---

#### Step 6  
Add the following paths to the **INCLUDE** environment variable:

```
(your VS path)\(xxxx)\BuildTools\VC\Tools\MSVC\(xx.xx.xxxxx)\include
(your WinKits path)\(xx)\Include\(xx.x.xxxxx.x)\ucrt
(your WinKits path)\(xx)\Include\(xx.x.xxxxx.x)\um
(your WinKits path)\(xx)\Include\(xx.x.xxxxx.x)\winrt
(your WinKits path)\(xx)\Include\(xx.x.xxxxx.x)\cppwinrt
(your WinKits path)\(xx)\Include\(xx.x.xxxxx.x)\shared
```

After configuration, restart your terminal and verify again:

```bash
cl
```

---

### 3️⃣ Install `multireflects`

Once all dependencies and the compiler are properly configured, run:

```bash
python setup.py
```

The package should now install successfully.

---

## 📖 Usage

See the file `instructions.pdf`.
