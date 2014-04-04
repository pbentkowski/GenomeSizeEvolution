#!/usr/bin/env python
import pylab as p
import matplotlib.cm as cm

# ---- importing and calculating data ---
data = p.genfromtxt("GeneralData.dat")
genome_size = p.genfromtxt("GenomeSizeData.dat")
age_data = p.genfromtxt("CellsAgeData.dat")
rel_genome_size = p.zeros((genome_size.shape[0], genome_size.shape[1]))
rel_age_data = p.zeros((age_data.shape[0], age_data.shape[1]))

# calculating stuff for genome size distribution plot
for i in range(genome_size.shape[0]):
        S = sum(genome_size[i, :])
        for j in range(genome_size.shape[1]):
            rel_genome_size[i, j] = genome_size[i, j] / S

# calculating stuff for age of cells distribution plot
for k in range(age_data.shape[0]):
        AgeSum = sum(age_data[k, :])
        for l in range(age_data.shape[1]):
            rel_age_data[k, l] = age_data[k, l] / AgeSum

frame_half_size = 300
mean_shannon = p.zeros((data.shape[0]-(2*frame_half_size), 3))
mean_shannon[:, 0] = data[frame_half_size:-frame_half_size][:, 0].copy()
k = frame_half_size

for j in range(mean_shannon.shape[0]):
        mean_shannon[j, 1] = data[k - frame_half_size:k
                                  + frame_half_size, 9].mean()
        mean_shannon[j, 2] = data[k - frame_half_size:k
                                  + frame_half_size, 9].std()
        k = k + 1

# ---- end of importing and calculating data  ---
FZ = 12
#---- first plot-------
p.figure(1, figsize=(12, 6))

p.subplot(311)
p.plot(data[:, 0], data[:, 1], '-')
p.axis([0, data[:, 0].max(), -1, 1])
p.ylabel('environmental conditions', fontsize=FZ)
p.grid(True)

p.subplot(312)
p.plot(data[:, 0], data[:, 4], '-')
p.axis([0, data[:, 0].max(), 0, data[:, 4].max()])
p.ylabel('number of cells', fontsize=FZ)
p.grid(True)

p.subplot(313)
p.plot(data[:, 0], data[:, 2], '-', linewidth=1)
p.fill_between(data[:, 0], data[:, 2] + data[:, 3], data[:, 2] - data[:, 3],
               color=(0.75, 0.75, 1.0, 0.75))
p.axis([0, data[:, 0].max(), 0, (data[:, 2] + data[:, 3]).max()])
p.xlabel('time (steps)', fontsize=FZ)
p.ylabel('mean number of genes', fontsize=FZ)
p.grid(True)
#---- end of first plot-------

#---- secont plot-------
p.figure(2, figsize=(12, 6))

p.subplot(311)
p.plot(data[:, 0], data[:, 5], '-')
p.axis([data[:, 0].min(),  data[:, 0].max(), 0,  data[:, 5].max()])
p.ylabel('resources inside cells')
p.grid(True)

p.subplot(312)
p.plot(data[:, 0], data[:, 6], '-')
p.axis([data[:, 0].min(),  data[:, 0].max(), 0,  data[:, 6].max()])
p.ylabel('resources in the env.')
p.grid(True)

p.subplot(313)
p.plot(data[:, 0], data[:, 5]+data[:, 6], '-')
p.axis([data[:, 0].min(),  data[:, 0].max(), (data[:, 5]+data[:, 6]).min(),
        (data[:, 5]+data[:, 6]).max()])
p.xlabel('time (steps)')
p.ylabel('diff. res.')
p.grid(True)
#---- end of second plot-------

#---- third plot-------
x = data[:, 0]
y = p.r_[1:genome_size.shape[1]+1]
p.figure(3, figsize=(12, 6))
p.contourf(x, y, rel_genome_size.transpose(), cmap=cm.jet)
p.axis([data[:, 0].min(),  data[:, 0].max(),  y.min(),  y.max()])
p.xlabel('time (steps)', fontsize=FZ)
p.ylabel('number of genes', fontsize=FZ)
p.grid(True)
cb = p.colorbar(format=r"%.2f")
cb.set_label('frequency of occurrence', fontsize=FZ)
#---- end of third plot-------


#---- fourth plot-------
p.figure(4, figsize=(10, 6))

p.subplot(311)
p.plot(data[:, 0], data[:, 11], '-')
p.axis([data[:, 0].min(),  data[:, 0].max(), 0.9 * data[:, 11].min(),
        1.1*data[:, 11].max()])
p.ylabel('number of dead')
p.grid(True)

p.subplot(312)
p.plot(data[:, 0], data[:, 10], '-')
p.axis([data[:, 0].min(),  data[:, 0].max(), 0.9*data[:, 10].min(),
        1.1*data[:, 10].max()])
p.ylabel('number of born')
p.grid(True)

p.subplot(313)
p.plot(data[:, 0], data[:, 9], '-')
p.plot(mean_shannon[:, 0], mean_shannon[:, 1], 'r-',
       label='stepping mean', linewidth=2)
p.axis([data[:, 0].min(),  data[:, 0].max(), 0.9 * data[:, 9].min(),
        1.1*data[:, 9].max()])
p.xlabel('time (steps)')
p.ylabel('Shannon index')
p.grid(True)
#---- end of fourth plot-------

#---- fifth plot-------
yy = p.r_[1:age_data.shape[1]+1]
p.figure(5, figsize=(12, 6))
p.contourf(x, yy, p.log(rel_age_data.transpose()), cmap=cm.jet)
#p.contourf(10.*x, yy, rel_age_data.transpose(), cmap=cm.jet)
p.axis([data[:, 0].min(),  data[:, 0].max(),  yy.min(),  yy.max()])
p.xlabel('time (steps)')
p.ylabel('Age of cells (steps)')
p.grid(True)
cbb = p.colorbar(format=r"%.2f")
cbb.set_label('Logarithm of frequency of occurrence')
#---- end of fifth plot-------

p.show()
print "done!"
