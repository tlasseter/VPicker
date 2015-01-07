package com.Sts.Actions.Wizards.Collaboration;

import com.Sts.IO.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsCollaborationParticipant
{
	public String username;
	public String name;
	public String ipAddress;
	public int minPort;
	public int maxPort;

	public StsCollaborationParticipant(String username, String name, String ipAddress, int minPort, int maxPort)
	{
		this.username = username;
		this.name = name;
		this.ipAddress = ipAddress;
		this.minPort = minPort;
		this.maxPort = maxPort;
	}
	static public StsCollaborationParticipant[] readFile(String directory, String filename)
	{
		String[] strings;

		StsFile file = StsFile.constructor(directory, filename);
		try
		{
			strings = file.readStringsFromFile();
			int nStrings = strings.length;
			StsCollaborationParticipant[] participants = new StsCollaborationParticipant[0];
			for (int n = 0; n < nStrings; )
			{
				String username = strings[n++];
				String name = strings[n++];
				String ipAddress = strings[n++];
				int minPort = Integer.parseInt(strings[n++]);
				int maxPort = Integer.parseInt(strings[n++]);

				participants = (StsCollaborationParticipant[])StsMath.arrayAddElement(participants, new StsCollaborationParticipant(username, name, ipAddress, minPort, maxPort));
			}
			return participants;
		}
		catch (Exception e)
		{
			StsException.systemError("Failed to read strings from" + directory + filename);
			return null;
		}
	}

	static StsCollaborationParticipant getParticipantWithUsername(StsCollaborationParticipant[] participants, String username)
	{
		for(int n = 0; n < participants.length; n++)
			if(participants[n].username.equals(username)) return participants[n];
		return null;
	}

	static StsCollaborationParticipant getParticipantWithName(StsCollaborationParticipant[] participants, String name)
	{
		for(int n = 0; n < participants.length; n++)
			if(participants[n].name.equals(name)) return participants[n];
		return null;
	}

	public String toString() { return name; }
}
