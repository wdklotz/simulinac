def deco(f):
    print('deco({}) called'.format(f.__name__))
    def wrapper(arg):
        print('in wrapper(...) before call to {}(...)'.format(f.__name__))
        f(arg.upper())
        print('in wrapper(...) after calling {}(...)'.format(f.__name__))
    return wrapper

@deco
def func(text):
    print('bare function func(...) called with argument "{}"'.format(text))

@deco
def greeting(message):
    print('greeting: "{}"'.format(message))

func('decorated function')
func('hello world')
greeting('happy coding')