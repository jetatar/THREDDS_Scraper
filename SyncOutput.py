#
#   Usage: python SyncOutput.py [source] [destination] [periodicity (in sec)]
#

import os, subprocess, shlex, datetime, time, sys


if len( sys.argv ) < 5:
	exit( 'Please provide four arguments - [Source] [Destination] [Periodicity]' )

source  = sys.argv[1]
target  = sys.argv[2]
mode    = sys.argv[3]
period  = int( sys.argv[4] )

if not source.endswith( '/' ):
	source += '/'
if not target.endswith( '/' ):
	target += '/'

mode = mode.lower()	# Lower case input


while True:

	cmd = 'rsync -avzi --safe-links'

	print( "++ %s:  Syncing all content." % datetime.datetime.now() )

	cmd += ' %s %s' % (source, target)
	print( "++ %s:  Executing: %s" % (datetime.datetime.now(), cmd) )
	command = shlex.split( cmd )
	error = subprocess.Popen( command, stderr = subprocess.PIPE ).communicate()

	if( int(period) != 0 ):
		print( "++ %s:  Sleeping for %s seconds." \
										% (datetime.datetime.now(), period) )
		time.sleep( int(period) )
