#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 19:02:36 2010

@author: Piotr Bentkowski
"""
import pylab as p
import gaussar

FontSize = 20

x = p.arange(-1.0, 1.0, 0.01)
genes = p.array([[1.0, 0.06, -0.7], [0.5, 0.3, 0.3], [0.12, 0.5, -0.1]])
y = p.zeros((x.shape[0], genes.shape[0]))
i = 0
for item in genes:
    y[:, i] = gaussar.gaussar(x, item)
    i += 1


p.figure(1, figsize=(10, 8))
p.plot(x, y, 'b-', linewidth=2)
p.xlabel('environmental conditions [ $x$ ]', fontsize=FontSize)
p.ylabel('efficiemcy of resource uptake [ $u(x,c,\sigma,A)$ ]',
         fontsize=FontSize)
p.text(0.25, 0.85, r"$u(x,c,\sigma,A)= A e^{\frac{-(x-c)^2}{2\sigma^2}} $",
       fontsize=26)
p.text(0.25, 0.75, r"$\sigma = \frac{A \alpha}{\sqrt{2\pi}} $", fontsize=26)
p.text(0.24, 0.25, r"$A $", fontsize=25)
p.text(0.4, 0.305, r"$\sigma $", fontsize=25)
p.text(0.24, 0.01, r"$c $", fontsize=25)
p.arrow(0.3, 0.0, 0.0, 0.465, head_width=0.015, head_length=0.03, color='red')
p.arrow(0.3, 0.5, 0.0, -0.466, head_width=0.015, head_length=0.03, color='red')
p.arrow(0.6, 0.3, -0.27, 0.0, head_width=0.015, head_length=0.03, color='red')
p.arrow(0.3, 0.3, 0.27, 0.0, head_width=0.015, head_length=0.03, color='red')
p.xticks(size=19)
p.yticks(size=19)
p.show()
