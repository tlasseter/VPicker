package com.Sts.Collaboration;

import java.io.*;
/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class TransactionMessage implements Serializable
{
	byte[] transactionBytes;
	String description;
	byte transactionType;
    public TransactionMessage(byte[] bytes, String description, byte transactionType)
    {
		transactionBytes = bytes;
		this.description = description;
		this.transactionType = transactionType;
    }
}
