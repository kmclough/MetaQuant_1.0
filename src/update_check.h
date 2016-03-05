/**
 *  update_check.h
 *  MetaQuant
 *
 *  Kevin McLoughlin
 *  Based on code from eXpress, created by Adam Roberts in 2011.
 *  Based on code from http://www.linuxhowtos.org/C_C++/socket.htm
 *  Modified by Adam Roberts on 1/18/11.
 *  Copyright 2014 Kevin McLoughlin, Adam Roberts. All rights reserved.
 */

// =-= TODO: Update URLs embedded in this code and reenable version checking in main.cpp,
// =-= if we ever have a dedicated web site for MetaQuant.

#include <signal.h>
#include <strings.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

int NUM_SEPS = 3;
int CONNECT_TIMEOUT = 5;

static int sTimeout = 0;

static void AlarmHandler(int sig)

{

	sTimeout = 1;

}

bool error(const char *msg)
{
	return false;
}

int parse_version_str(char* version_str)
{
	int version_int = 0;
	char* token = strtok(version_str,".");
    for(int i = 0; i < NUM_SEPS; ++i)
    {
      version_int += atoi(token)*(int)pow(100.,NUM_SEPS-i);
	}
	return version_int;
}

bool get_current_version(char* curr_version)
{
    int sockfd, portno, n;
    struct sockaddr_in serv_addr;
    struct hostent *server;
	
    portno = 80;
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0)

        return error("ERROR opening socket");
	
    server = gethostbyname("bio.math.berkeley.edu");
    if (server == NULL)

        return error("ERROR, no such host");

    bzero((char *) &serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    bcopy((char *)server->h_addr,

		  (char *)&serv_addr.sin_addr.s_addr,
		  server->h_length);
    serv_addr.sin_port = htons(portno);

	signal(SIGALRM, AlarmHandler);

	sTimeout = 0;

	alarm(CONNECT_TIMEOUT);

	
	int ret;
	ret = connect(sockfd, (struct sockaddr*)&serv_addr, sizeof(serv_addr));
	if (ret < 0 || sTimeout)
	{
		return error("ERROR connecting");
	}
	
	char buffer[1024];
	strcpy(buffer, "GET /eXpress/curr_xprs_version HTTP/1.1\nHost: bio.math.berkeley.edu\n\n");
	n = (int)write(sockfd,buffer,1024);
	
    if (n < 0)

		return error("ERROR writing to socket");
	bzero(curr_version, sizeof(curr_version));
    n = (int)read(sockfd,buffer,1024);
    if (n < 0)

		return error("ERROR reading from socket");

	char* token;
	token = strtok(buffer, "$");
	token = strtok(NULL, "$");
	if (token==NULL)
		return error("ERROR parsing response");
	
	strcpy(curr_version, token);
		
	return true;
}

void check_version(const char* this_version)
{
	char curr_version[256];
    memset(curr_version, 0, sizeof(curr_version));
	if (get_current_version(curr_version))
	{
		if (strcmp(curr_version, this_version)==0) {
			fprintf(stderr, "You are using MetaQuant v%s, which is the most recent release.\n\n", PACKAGE_VERSION);
		} else {
			fprintf(stderr, "WARNING: Your version of MetaQuant is not up-to-date. It is recommended that you upgrade to v%s to benefit from the most recent features and bug fixes.\n\n", curr_version);
    }
	} else {
		fprintf(stderr, "WARNING: Could not connect to update server to verify current version. Please check at the MetaQuant website.\n\n");
	}
}
