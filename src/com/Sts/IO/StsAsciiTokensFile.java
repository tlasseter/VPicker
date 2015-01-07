package com.Sts.IO;

import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jan 24, 2008
 * Time: 9:47:34 PM
 * To change this template use File | Settings | File Templates.
 */

/* This ascii file consists of a set of lines, each line containing a defined number of tokens */
public class StsAsciiTokensFile extends StsAsciiFile
{
    /** if nTokensPerLine is not initialized to a postive number, there are assumed to be a variable and/or unknown number of tokens per line */
    int nTokensPerLine = 0;
    public String[][] tokenTable = new String[nPropertiesInc][];
    public int nTableLines = 0;
    String[] delimiters;

    static final String spaceDelimiter = " ";
    
    static int nPropertiesInc = 10;

    public StsAsciiTokensFile(String directory, String filename, int nTokensPerLine, String delimiter)
    {
        file = StsFile.constructor(directory, filename);
        initialize(file);
        this.nTokensPerLine = nTokensPerLine;
        if (delimiter == null) delimiter = " ";
        delimiters = new String[] { delimiter };
        readFile();
    }

    public StsAsciiTokensFile(String pathname, int nTokensPerLine, String delimiter)
    {
        file = StsFile.constructor(pathname);
        initialize(file);
        this.nTokensPerLine = nTokensPerLine;
        delimiters = new String[] { delimiter };
        readFile();
    }

    public StsAsciiTokensFile(String directory, String filename, int nTokensPerLine)
    {
        this(directory, filename, nTokensPerLine, spaceDelimiter);
    }

    public StsAsciiTokensFile(String directory, String filename)
    {
        this(directory, filename, 0, null);
    }

    public StsAsciiTokensFile(String directory, String filename, String delimiter)
    {
        this(directory, filename, 0, delimiter);
    }

    static public String[][] getStringTable(String directory, String filename, String delimiter)
    {
        try
        {
            StsAsciiTokensFile file = new StsAsciiTokensFile(directory, filename, delimiter);
            return file.tokenTable;
        }
        catch(Exception e)
        {
            return null;
        }
    }

    static public String[][] getStringTable(String pathname, String delimiter)
    {
        try
        {
            StsAsciiTokensFile file = new StsAsciiTokensFile(pathname, delimiter);
            return file.tokenTable;
        }
        catch(Exception e)
        {
            return null;
        }
    }

    /** read the tokens from the file if it exists; if it doesn't exist, create an empty file */
    static String[][] readAsciiFileTokens(String directory, String filename, int nTokensPerLine)
    {
        try
        {
            StsAsciiTokensFile tokensFile = new StsAsciiTokensFile(directory, filename, nTokensPerLine);
            return tokensFile.readFile();
        }
        catch(Exception e)
        {
            StsException.outputException("StsAsciiTokensFile.getAsciiFileTokens() failed for: " +
                    directory + File.separator + filename, e, StsException.WARNING);
            return null;
        }
    }

    public String[][] readFile()
    {
        String[] tokens;
        try
        {
            if(!openRead()) return null;
            while (true)
            {
                if(delimiters != null)
                    tokens = getTokens(delimiters);
                else
                    tokens = getTokens();
                if(tokens == null) break;
                addLineTokens(tokens);
            }
            tokenTable = (String[][])StsMath.trimArray(tokenTable, nTableLines);
            return tokenTable;
        }
        catch (Exception e)
        {
            new StsMessage(null, StsMessage.WARNING, "File read error for " +
                           getFilename() + ": " + e.getMessage());
            return null;
        }
        finally
        {
            close();
        }
    }

    public boolean writeFile()
    {
        if(!openWriteWithErrorMessage()) return false;
        try
        {
            for(int n = 0; n < nTableLines; n++)
            {
                String[] tokens = tokenTable[n];
                for(int i = 0; i < tokens.length; i++)
                    writer.write(tokens[i] + spaceDelimiter);
                writer.write("\n");
            }
            close();
            return true;
        }
        catch(Exception e)
        {
            StsException.systemError(this, "writeFile", "couldn't write file: " + getFilename());
            return false;
        }
        finally
        {
            close();
        }
    }

    public boolean addLineTokens(String[] tokens)
    {
        int nTokens = tokens.length;
        if(nTokensPerLine != 0 && nTokens != nTokensPerLine)
        {
            StsException.systemError(this, "addLineTokens", "failed to find " + nTokensPerLine + " tokens on line " + getLine());
            return false;
        }
        int length = tokenTable.length;
        if(tokenTable == null || tokenTable.length == nTableLines)
        {
            String[][] newProperties = new String[length + nPropertiesInc][];
            System.arraycopy(tokenTable, 0, newProperties, 0, length);
            tokenTable = newProperties;
        }
        tokenTable[nTableLines++] = tokens;
        return true;
    }

    public String[] getMatchToToken(String match)
    {
        if(nTableLines == 0) return null;
        for(int n = 0; n < nTableLines; n++)
        {
            if(match.equals(tokenTable[n][0]))
                return tokenTable[n];
        }
        return null;
    }
    
    public String[][] getMatchToTokens(String match)
    {
        if(nTableLines == 0) return null;
        String[][] matchTokens = new String[tokenTable.length][];
        int cnt = 0;
        for(int n = 0; n < nTableLines; n++)
        {
            if(match.equals(tokenTable[n][0]))
                matchTokens[cnt++] = tokenTable[n];
        }
        if(cnt > 0)
        	return (String[][])StsMath.trimArray(matchTokens, cnt);
        else
        	return null;
    }
    
    public void debugPrint()
    {
        for(int n = 0; n < nTableLines; n++)
        {
            System.out.print(n + ":");
            for(int i = 0; i < nTokensPerLine; i++)
                System.out.print(spaceDelimiter + tokenTable[n][i]);
            System.out.println();
        }
    }

    static public void main(String[] args)
    {
        try
        {
            StsModel.constructor();
            // String directory = System.getProperty("user.home") + File.separator + "S2SCache" + File.separator;
            // String filename = "directoryAliases.txt";
            String directory = "c:\\data\\baker\\DTS\\DTS data (corrected)";
            String filename = "P2S2 corrected DTS.csv";
            StsAsciiTokensFile tokenStringsFile = new StsAsciiTokensFile(directory, filename, ", ");
        /*
            for(int n = 0; n < 5; n++)
            {
                String[] newStringTokens = new String[] { "initial" + n, "initial" + n + "A"}; 
                tokenStringsFile.addLineTokens(newStringTokens);
            }
            tokenStringsFile.writeFile();
        */
         //   tokenStringsFile.readFile();
            tokenStringsFile.debugPrint();
            /*
            String[] newStringTokens = new String[] { "added1", "added1A"};
            tokenStringsFile.addLineTokens(newStringTokens);
            tokenStringsFile.writeFile();
        */
        }
        catch(Exception e)
        {
            StsException.outputException("StsAsciiTokensFile.main() failed.", e, StsException.WARNING);
        }
    }

    /**
     * sets first instance in table that matches "key" to "value"
     * if key is not in table yet, it will be added.
     * @param key
     * @param value
     */
    public boolean setToken(String key, String value)
    {
        for (int i=0; i< tokenTable.length; i++)
        {
            if (tokenTable[i] != null && tokenTable[i][0].equals(key))
            {
                tokenTable[i][1] = value;
                return true;
            }
        }
        return addLineTokens(new String[]{key, value});  //if here, key doesn't exist yet so token needs to be added
    }
}
