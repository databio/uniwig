import itersum
import time

# Here's a python generator function. Python generators are functions that return
# a sequence of results. You can create a generator by just writing a function,
# but using the `yield` keyword instead of `return`. This acts like a 'return' keyword,
# but the function remembers where it was so that the next time it is called, it
# picks up where it left off.
# More info on generators: https://wiki.python.org/moin/Generators

# Calling this generator function will yield a *generator object*, which behaves
# like an iterator, with the key advantage that it enumerates items on the fly,
# rather than storing them all in memory.

# Generator function
def countdown(num):
    print('Starting')
    if num == 5000:
    	time.sleep(5)
    while num > 0:
        yield num
        num -= 1


i = countdown(1500)

# Create an iterator by calling the countdown function

itersum.itersum(countdown(1500))
itersum.itersum(countdown(4999))
itersum.itersum(countdown(5001))

