package com.Sts.Utilities;

import com.Sts.Interfaces.*;

import javax.swing.*;
import java.io.*;
import java.util.*;

/**
 * Title:
 * Description:
 * Copyright:    Copyright (c) 2001
 * Company:
 * @author
 * @version 1.0
 */

public class StsRandomAccessFileListModel extends AbstractListModel implements StsListModel
{
    int type;
	RandomAccessFile file;
	int maxLinesPerBlock;
	ArrayList blocks = new ArrayList(10);
	Block lastBlock;
	BlocksFIFO blocksFIFO = new BlocksFIFO(3);
	long endPos = -1;
	int endLine = -1;
	String longestString = "";

    public StsRandomAccessFileListModel(RandomAccessFile file, int maxLinesPerBlock, int type)
    {
 		String string;

		try
		{
            this.type = type;
			this.file = file;
			this.maxLinesPerBlock = maxLinesPerBlock;
			file.seek(0);
			while(readLine());
			endLine = lastBlock.endLine;
			endPos = file.getFilePointer()-1;
			clearBlockLines();
		}
		catch(Exception e)
		{
			System.out.println("BlockDocument failed.\n");
			e.printStackTrace();
		}
    }

	private void addBlock(long startPos, int startLine) throws IOException
	{
		lastBlock = new Block(startPos, startLine);
		blocks.add(lastBlock);
	}

	private void checkBlockIsFull()
	{
		long startPos = 0;
		int startLine = 0;

		try
		{
			if(lastBlock == null)
			{
				startPos = file.getFilePointer();
			    addBlock(startPos, startLine);
			}
			else if(!lastBlock.isFull())
				return;
			else
			{
			    startPos = lastBlock.endPos+1;
			    startLine = lastBlock.endLine+1;
			    addBlock(startPos, startLine);
			}
		}
		catch(Exception e)
		{
			System.out.println("StsRandomAccessFileListModel.checkBlockIsFull() failed.");
			e.printStackTrace();
		}
	}

	// reads lines from file into blocks
	public boolean readLine()
	{
		try
		{
			checkBlockIsFull();
			if(!lastBlock.readLine()) return false;
			endLine = lastBlock.endLine;
			fireIntervalAdded(this, endLine, endLine);
			return true;
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return false;
		}
	}

	public void addLine(String line)
	{
		writeLine(line);
	}

	public void writeLine(String line)
	{
		try
		{
			checkBlockIsFull();
			lastBlock.writeLine(line); // + lineEndString
			endLine = lastBlock.endLine;
			fireIntervalAdded(this, endLine, endLine);
		}
		catch(Exception e)
		{
//			System.out.println("StsRandomAccessFileListModel.writeLine() failed.\n" + " Line: " + line);
		}
	}

	public String getStringAtLine(int nLine)
	{
		try
		{
			Iterator iter = blocks.iterator();
			while(iter.hasNext())
			{
				Block block = (Block)iter.next();
				String string = block.getStringAtLine(nLine);
				if(string != null) return string;
			}
			return null;
		}
		catch(Exception e)
		{
			System.out.println("StsRandomAccessFileListModel.getStringAtLine() failed.");
			e.printStackTrace();
			return null;
		}
	}

	private void clearBlockLines()
	{
		Iterator iter = blocks.iterator();
		while(iter.hasNext())
		{
			Block block = (Block)iter.next();
			block.lines = null;
		}
	}

	public String getLongestString() { return longestString; }

	// Implementation of ListModel interface
    public int getSize()
    {
		return endLine+1;
    }
    public Object getElementAt(int index)
    {
		return getStringAtLine(index);
    }

	public void closeFile()
	{
		try
		{
		    if(file != null) file.close();
		}
		catch(Exception e)
		{
//			StsException.outputException("StsRandomAccessFileListModel.closeFile() failed.",
//				e, StsException.WARNING);
		}
	}

	class Block
	{
		int nLines = 0;
		long startPos, endPos;
		int startLine, endLine;
		String[] lines;

		Block(long startPos, int startLine) throws IOException
		{
			this.startPos = startPos;
			this.startLine = startLine;
			endPos = startPos;
			endLine = startLine-1;
			lines = new String[maxLinesPerBlock];
		}

		void writeLine(String string) throws IOException
		{
		    if(lines == null) restoreLines();
			lines[nLines++] = string;
			checkForLongestLine(string);
			file.seek(endPos+1);
			file.writeBytes(string + "\n");
			endPos = file.getFilePointer()-1;
			endLine++;
		}

		boolean readLine() throws IOException
		{
		    if(lines == null) return false; // should not happen
		    String string = file.readLine();
			if(string == null) return false;
			lines[nLines++] = string;
			checkForLongestLine(string);
			endPos = file.getFilePointer()-1;
			endLine++;
			return true;
		}

		void checkForLongestLine(String string)
		{
			if(longestString.length() < string.length()) longestString = string;
		}

		boolean isFull()
		{
		    return nLines == maxLinesPerBlock;
		}

		String getStringAtLine(int nLine)
		{
			if( nLine < startLine || nLine - startLine >= maxLinesPerBlock) return null;
			if(lines == null) restoreLines();
			return lines[nLine-startLine];
		}

		void restoreLines()
		{
			try
			{
				file.seek(startPos);
				lines = new String[maxLinesPerBlock];
				for(int n = 0; n < nLines; n++)
				{
					String line = file.readLine();
					lines[n] = line;
				}
				blocksFIFO.add(this);
			}
			catch(Exception e)
			{
				System.out.println("StsRandomAccessFileListModel.restoreLines() failed.");
				e.printStackTrace();
			}
		}
	}

	// maintain an array of blocks whose lines arrays are non-null.
	// All other blocks have their lines arrays nulled out.
	// If we need the lines for a nulled block, remove the top one (the oldest),
	// null it out, and add the new one at the bottom.

	class BlocksFIFO
	{
		Block[] blocks;
		int maxSize;
		int nBlocks = 0;

	    BlocksFIFO(int maxSize)
		{
			this.maxSize = maxSize;
			blocks = new Block[maxSize];
		}

		void add(Block block)
		{
			if(nBlocks == maxSize)
			{
				blocks[0].lines = null;  // remove lines array from this block
				for(int n = 1; n < maxSize; n++)
				    blocks[n-1] = blocks[n];
				nBlocks--;
			}
			blocks[nBlocks++] = block;
		}
	}
}
