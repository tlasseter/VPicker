package com.Sts.DB;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import javax.jnlp.*;
import java.io.*;

public class StsRandomAccessFile extends RandomAccessFile implements JNLPRandomAccessFile
{
	String pathname;

	public StsRandomAccessFile(File file, String mode) throws FileNotFoundException
	{
		super(file, mode);
		pathname = file.getPath();
	}

	public StsRandomAccessFile(String pathname, String mode) throws FileNotFoundException
	{
		super(pathname, mode);
		this.pathname = pathname;
	}

	public String getPathname() { return pathname; }
}
