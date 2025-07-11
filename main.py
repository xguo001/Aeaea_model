# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    start = 0
    alive = True

    while alive:
        end = start
        start += 1
        if start == 10:
            alive = False

        print ("start", start)
        print ("end", end)


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
