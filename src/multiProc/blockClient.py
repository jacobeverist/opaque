# This is the blocking Get Poetry Now! client.

import datetime, optparse, socket

def parse_args():
    usage = """usage: %prog [options] [hostname]:port ...

This is the Get Poetry Now! client, blocking edition.
Run it like this:

  python get-poetry.py port1 port2 port3 ...

If you are in the base directory of the twisted-intro package,
you could run it like this:

  python blocking-client/get-poetry.py 10001 10002 10003

to grab poetry from servers on ports 10001, 10002, and 10003.

Of course, there need to be servers listening on those ports
for that to work.
"""

    parser = optparse.OptionParser(usage)

    _, addresses = parser.parse_args()

    if not addresses:
        print parser.format_help()
        parser.exit()

    def parse_address(addr):
        if ':' not in addr:
            host = '127.0.0.1'
            port = addr
        else:
            host, port = addr.split(':', 1)

        if not port.isdigit():
            parser.error('Ports must be integers.')

        return host, int(port)

    return map(parse_address, addresses)

def send_problem(address, data):
    """Download a piece of poetry from the given address."""

    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect(address)

    sock.sendall(data + '\n')

    results = ''

    while True:

        # This is the 'blocking' call in this synchronous program.
        # The recv() method will block for an indeterminate period
        # of time waiting for bytes to be received from the server.

        bytes = sock.recv(1024)

        if not bytes:
            sock.close()
            break

        results += bytes

    return results

def format_address(address):
    host, port = address
    return '%s:%s' % (host or '127.0.0.1', port)


def main():
	f = open("shapeInitData.txt")
	data = eval(f.read().rstrip())
	f.close()
	results = eval(send_problem(('127.0.0.1', 10000), repr(data)))
	print len(results), "elements"
	print results[0]

	#addresses = parse_args()

	#for i, address in enumerate(addresses):
	#    addr_fmt = format_address(address)

	#    print 'Problem solver: %s' % (addr_fmt)


if __name__ == '__main__':
    main()
