#!/usr/bin/env python3
#
# Expectation Maximization for Gaussian Mixture Models.
# Copyright (C) 2012-2022 Juan Valor
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import sys,random,math

# Write to disk the cache data and update the console message.
def write_cache(actual,samples,cache,file):
	file.write("\n".join(cache)+"\n")
	sys.stdout.write("\rGenerated sample "+str(actual)+" of "+str(samples)+".")
	return 0,[]

# Generate random Gaussian numbers based on the distributions.
def generate(samples,dimensions,distributions,file):
	count,cache=write_cache(0,samples,[str(dimensions)+" "+str(samples)],file)
	for i in range(1,samples+1):
		distribution,numbers=random.choice(distributions),[]
		for j in range(dimensions):
			m,v=distribution[j]
			numbers.append(("%.3f"%random.gauss(m,v)))
		if count>=100000: # Use a cache to speedup the write on disk.
			count,cache=write_cache(i,samples,cache,file)
		count+=dimensions
		cache.append(" ".join(numbers))
	count,cache=write_cache(samples,samples,cache,file)
	sys.stdout.write("\n")

# Generate some random Gaussian distributions (mean and covariance).
def generate_distributions(samples,dimensions):
	numdist,distributions=int(math.sqrt(samples/2.0))+1,[]
	for i in range(numdist):
		distdim=[]
		for j in range(dimensions):
			m,v=random.uniform(-10000,10000),random.uniform(10,100)
			distdim.append((m,v))
		distributions.append(distdim)
	return distributions

# The main execution of the random data generator.
def main():
	if len(sys.argv)<4:
		sys.stderr.write("Usage: "+sys.argv[0]+" <dimensions> <samples> <file.txt>\n")
	else:
		dimensions,samples=int(sys.argv[1]),int(sys.argv[2])
		file=open(str(sys.argv[3]),"w")
		distributions=generate_distributions(samples,dimensions)
		generate(samples,dimensions,distributions,file)
		file.close()

# If is not imported as a library, execute the main function.
if __name__=="__main__": main()

