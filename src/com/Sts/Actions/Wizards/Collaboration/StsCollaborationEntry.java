package com.Sts.Actions.Wizards.Collaboration;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

import com.Sts.IO.*;
import com.Sts.Utilities.*;

public class StsCollaborationEntry
{
	public String name;
	public String description;
	public StsCollaborationParticipant[] participants;

	public StsCollaborationEntry(String name, String description, StsCollaborationParticipant[] participants)
	{
		this.name = name;
		this.description = description;
		this.participants = participants;
	}

	static public StsCollaborationEntry[] readFile(String directory, String filename, StsCollaborationParticipant[] allParticipants)
	{
		String[] strings;

		StsFile file = StsFile.constructor(directory, filename);
		try
		{
			strings = file.readStringsFromFile();
			int nStrings = strings.length;
			StsCollaborationEntry[] entries = new StsCollaborationEntry[0];
			int nEntry = 0;
			for (int n = 0; n < nStrings; )
			{
				String name = strings[n++];
				String description = strings[n++];
				StsCollaborationParticipant[] participants = new StsCollaborationParticipant[0];
				while(n < nStrings)
				{
					String string = strings[n++];
					if(string.equals("")) break;
					StsCollaborationParticipant participant = StsCollaborationParticipant.getParticipantWithName(allParticipants, string);
					participants = (StsCollaborationParticipant[])StsMath.arrayAddElement(participants, participant);
				}
				entries = (StsCollaborationEntry[])StsMath.arrayAddElement(entries, new StsCollaborationEntry(name, description, participants));
			}
			return entries;
		}
		catch (Exception e)
		{
			StsException.systemError("Failed to read strings from" + directory + filename);
			return null;
		}
	}

	public String toString() { return name; }

	public StsCollaborationParticipant getLeader() { return participants[0]; }

	public StsCollaborationParticipant getParticipantWithUsername(String username)
	{
		return StsCollaborationParticipant.getParticipantWithUsername(participants, username);
	}

	public boolean isLeader(StsCollaborationParticipant participant) { return participant == participants[0]; }

	public boolean equals(StsCollaborationEntry otherEntry)
	{
		if(!name.equals(otherEntry.name)) return false;
		if(participants == null && otherEntry.participants == null) return false;
		if(participants.length != otherEntry.participants.length) return false;
		for(int n = 0; n < participants.length; n++)
			if(!participants[n].name.equals(otherEntry.participants[n].name)) return false;
		return true;
	}
}
