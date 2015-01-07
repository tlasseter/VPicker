

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.IO;

import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;

public class StsAsciiFile
{
    // instance fields
    StsAbstractFile file = null;
    String filename = null;
    BufferedReader reader = null;
    BufferedWriter writer = null;
    String line = null;
    int nLines = 0;
    long nBytes = 0;
    InputStream is = null;

    static public boolean debug = false;

    public StsAsciiFile()
    {
    }

    public StsAsciiFile(StsAbstractFile file)
    {
        initialize(file);
    }

    public void initialize(StsAbstractFile file)
    {
        this.file = file;
        filename = file.getFilename();
    }

    /** If file doesn't exist, don't throw exception, simply return.
     *  Call this method when file existence is not necessary, but simply indicates info doesn't exist
     *  or is not otherwise available.
     * @return
     * @throws IOException
     */
    public boolean openRead()
    {
        try
        {
            if(reader != null) reader.close();
			if(!file.isAFile()) return false;
            is = file.getInputStream();
            if(is == null) return false;
            reader = new BufferedReader(new InputStreamReader(is));
            return true;
        }
        catch(Exception e)
        {
            return false;
        }
    }

    public boolean openReadWithErrorMessage()
    {
        try
        {
            if(reader != null) reader.close();
            is = file.getInputStream();
            reader = new BufferedReader(new InputStreamReader(is));
            return true;
        }
        catch(Exception e)
        {
            new StsMessage(null, StsMessage.WARNING, "Failed to open file " + file.getPathname());
            return false;
        }
    }

    public boolean openWrite() throws IOException
    {
        if(writer != null) writer.close();
        OutputStream os = file.getOutputStream(false);
        if(os == null) return false;
        writer = new BufferedWriter(new OutputStreamWriter(os));
        return true;
    }

    public boolean openWriteWithErrorMessage()
    {
        try
        {
            if(writer != null) writer.close();
            OutputStream os = file.getOutputStream(false);
            if(os == null) return false;
            writer = new BufferedWriter(new OutputStreamWriter(os));
            return true;
        }
        catch(Exception e)
        {
            new StsMessage(null, StsMessage.WARNING, "File not found " + file.getPathname());
            return false;
        }
    }

    public boolean openWriteAppend()
    {
        try
        {
            if(writer != null) writer.close();
            OutputStream os = file.getOutputStream(true);
            writer = new BufferedWriter(new OutputStreamWriter(os));
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsAsciiFile.openWrite() failed." +
                "Can't write: " + filename, e, StsException.WARNING);
            return false;
        }
    }
    public boolean openReadWrite() throws IOException
    {
        return openRead() && openWrite();
    }

    /** close this binary file */
    public boolean close()
    {
        try
        {
            if (writer != null)
            {
                writer.flush();
                writer.close();
                writer = null;
            }
            if (reader != null)
            {
                reader.close();
                reader = null;
            }
			return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsAsciiFile.close() failed."
                    + "Unable to close file " + filename, e, StsException.WARNING);
			return false;
        }
    }

    public String readLine() throws IOException
    {
        line = reader.readLine();
        if(debug) System.out.println("line " + nLines + ": " + line);
        if(line == null) return null;
        nLines++;
        nBytes += line.length();
        return line;
    }

    public String readNextNonBlankLine() throws IOException
    {
        line = reader.readLine();
        if(debug) System.out.println("line " + nLines + ": " + line);
        if(line == null) return null;
        nLines++;
        if(line.length() == 0)
            return readNextNonBlankLine();
        else
        {
            nBytes += line.length();
            return line;
        }
    }

    public String[] getTokens(String delimiters) throws IOException
    {
        line = reader.readLine();
        if(line == null) return null;
        line.trim();
        int nTokens = 0;
        StringTokenizer stringTokenizer = null;

        stringTokenizer = new StringTokenizer(line, delimiters);
        nTokens = stringTokenizer.countTokens();

        if(nTokens == 0) return null;

        nLines++;
        nBytes += line.length();
        String[] tokens = new String[nTokens];
        int n = 0;
        while(stringTokenizer.hasMoreTokens())
            tokens[n++] = stringTokenizer.nextToken();
        return tokens;
    }

    public String[] getTokens(String[] delimiters) throws IOException
    {
        line = reader.readLine();
        if(line == null) return null;
        line.trim();
        int nTokens = 0;
        StringTokenizer stringTokenizer = null;       
        for(int i=0; i< delimiters.length; i++)
        {
        	if(delimiters[i].equals(" "))
        		stringTokenizer = new StringTokenizer(line);    // tab, space, cr, lf
        	else
        		stringTokenizer = new StringTokenizer(line, delimiters[i]);
        	nTokens = stringTokenizer.countTokens();
        	if(nTokens > 1) break;
        }
        if(nTokens == 0) return null;

        nLines++;
        nBytes += line.length();
        String[] tokens = new String[nTokens];
        int n = 0;
        while(stringTokenizer.hasMoreTokens())
            tokens[n++] = stringTokenizer.nextToken();
        return tokens;
    }
    
    public String[] getTokens() throws IOException
    {
        line = reader.readLine();
        if(line == null) return null;
        line.trim();

        StringTokenizer stringTokenizer = new StringTokenizer(line);
        int nTokens = stringTokenizer.countTokens();
        if(nTokens == 0) return null;

        nLines++;
        nBytes += line.length();
        String[] tokens = new String[nTokens];
        int n = 0;
        while(stringTokenizer.hasMoreTokens())
            tokens[n++] = stringTokenizer.nextToken();
        return tokens;
    }

    public Iterator getTokenIterator()
    {
        try
        {
            return new TokenIterator();
        }
        catch(Exception e)
        {
            StsException.outputException("StsAsciiFile.getTokenIterator() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    class TokenIterator implements Iterator
    {
        StringTokenizer stringTokenizer;
        Object nextToken = null;

        TokenIterator() throws IOException
        {
            line = reader.readLine();
            nLines++;
            if(line == null) return;
            line = line.trim();
            stringTokenizer = new StringTokenizer(line);
        }

        public boolean hasNext()
        {
            try
            {
                if(stringTokenizer.hasMoreElements())
                {
                    nextToken = stringTokenizer.nextElement();
                    return true;
                }
                line = reader.readLine();
                nLines++;
                if(line == null) return false;
                stringTokenizer = new StringTokenizer(line);
                if(!stringTokenizer.hasMoreElements()) return false;
                nextToken = stringTokenizer.nextElement();
                return true;
            }
            catch(Exception e)
            {
                StsException.outputException("StsAsciiFile.TokenIterator.hasNext() failed.",
                    e, StsException.WARNING);
                return false;
            }
        }

        public Object next()
        {
            return nextToken;
        }

        public void remove()
        {
        }
    }

    public String getFilename() { return filename; }
    public String getLine() { return line; }
    public int getNLines() { return nLines; }
    public long length() { return file.length(); }
    public long getNBytes() { return nBytes; }

    public void writeString(String string) throws IOException
    {
        writer.write(string + " ");
    }

    public void writeLine(String string) throws IOException
    {
        writer.write(string);
        writer.newLine();
    }

    public void endLine() throws IOException
    {
        writer.newLine();
    }

    static public void main(String[] args)
    {
        StsAsciiFile asciiFile = null;
        try
        {
            String pathname = "c:\\asciiWriteTest.txt";
            StsFile stsFile = StsFile.constructor(pathname);
            asciiFile = new StsAsciiFile(stsFile);
            asciiFile.openWrite();
            asciiFile.writeLine("WELLNAME");
            asciiFile.writeLine("ORIGIN XY");
        }
        catch(Exception e)
        {
            StsException.outputWarningException(StsAsciiFile.class, "main", e);
        }
        finally
        {
            if(asciiFile != null) asciiFile.close();
        }
    }
}