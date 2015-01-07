package com.Sts.Collaboration;

import java.io.*;

class RequestDBMessage implements Serializable
{
	Participant p;
	RequestDBMessage(Participant p)
	{
		this.p = p;
	}
}
