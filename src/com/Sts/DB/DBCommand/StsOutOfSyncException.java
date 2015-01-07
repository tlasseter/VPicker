package com.Sts.DB.DBCommand;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsOutOfSyncException extends RuntimeException
{
	public StsOutOfSyncException()
	{
		super();
	}

	public StsOutOfSyncException(String message)
	{
		super(message);
	}

	public StsOutOfSyncException(Throwable cause)
	{
		super(cause);
	}

	public StsOutOfSyncException(String message, Throwable cause)
	{
		super(message, cause);
	}
}
