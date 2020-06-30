build-docker-image:
	sudo docker image build -t vadr:1.1 .

test-docker-image:
	sudo docker run vadr:1.1 '/bin/bash' '-c' '$VADRSCRIPTSDIR/testfiles/do-install-tests-local.sh'

