# Virtual envs on OSX

http://sourabhbajaj.com/mac-setup/Python/virtualenv.html

```bash
$ virtualenv venv --distribute --system-site-packages
```

These commands create a venv subdirectory in your project where everything is installed. You need to activate it first 
though (in every terminal where you are working on your project):

```bash
$ source venv/bin/activate
```

You should see a (venv) appear at the beginning of your terminal prompt indicating that you are working inside the 
virtualenv. Now when you install something:

```bash
$ pip install <package>
```

It will get installed in the venv folder, and not conflict with other projects.

To leave the virtual environment use.

```bash
$ deactivate
```
