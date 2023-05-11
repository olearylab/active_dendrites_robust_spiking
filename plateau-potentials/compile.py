"""
authors:
Salvador Dura <salvadordura@gmail.com>
Joe Graham <joe.w.graham@gmail.com>
"""
import os
import sys
from inspect import getsourcefile
import shutil
import subprocess

#h.load_file('stdrun.hoc')
# Get the current path
eeedir = os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))
# Get into the folder that containing all the mod files
moddir = os.path.abspath(os.path.join(eeedir, "mod"))

def compile(compiler="nrnivmodl", recompile=True):

	"""Removes any existing x86_64/i386 dirs, compiles using mkmod, and creates a symlink from sim dir to x86_64"""

	def rmcompdir(comppath=moddir):
		"""Removes x86_64 compiled mod directory from comppath dir."""
		path = os.path.join(comppath, "x86_64")
		if os.path.isdir(path):
			print("Removing directory: " + path)
			if (os.path.realpath(path) != path):
				os.remove(path)
			else:
				shutil.rmtree(path, ignore_errors=True)

	def compilemod(compiler=compiler):
		"""Compiles mod files using command line.
		Default compiler is 'nrnivmodl' (NeurosimLab-specific script). """
		print("Compiling mod file for EEE_Penny project using '" + compiler + "'...")
		compile_output = subprocess.call(compiler, shell=True)


	def linkmod(mydir):
		"""Symlinks the NEURON compiler output folder into the given directory."""
		if not os.path.isdir(os.path.join(moddir, "x86_64")):
			print("Compiled folder not found for symlinking...")
		else:
		#	if os.path.isdir(os.path.join(moddir, "x86_64")):
			if not os.path.isdir(os.path.join(mydir, "x86_64")):
				os.symlink(os.path.join(moddir, "x86_64"), os.path.join(eeedir, "x86_64"))
	#			if not os.path.isdir(os.path.join(batchdir, "x86_64")):
#					os.symlink(os.path.join(moddir, "x86_64"), os.path.join(batchdir, "x86_64"))


	curdir = os.getcwd()
	if curdir != eeedir:
		print("Moving to " + eeedir + " in order to compile mod files.")
		os.chdir(eeedir)

	rmcompdir(comppath=moddir)
	rmcompdir(comppath=eeedir)

	os.chdir(moddir)
	compilemod()
	linkmod(eeedir)
	os.chdir(eeedir)

	print("Compiling completed.")

###########################################

if __name__ == "__main__":
	print("Compiling the mod files")
	#make_output_dirs(all_batches)
	compile()
	print("Finished setting up EEE project.")
