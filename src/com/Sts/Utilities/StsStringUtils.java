package com.Sts.Utilities;

import java.util.StringTokenizer;


/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

public class StsStringUtils
{
    static public char[] invalidChars = new char[] {' ',':',';','"','\'',',', '[', ']', '{', '}', '=', '*', '@', '?','/', '\\', '$', '%', '&', '*'};
    // ')', '(', Removed.  Should not be included as invalid.
    static public String cleanString(String name)
    {
        return cleanString(name, invalidChars);
    }

    static public String cleanString(String name, char[] invalidChars)
    {
        for(int i=0; i<invalidChars.length; i++)
            name = name.replace(invalidChars[i],'_');

        return name;
    }

    static public String[] stringListAdd(String[] strings, String string)
    {
        return (String[])StsMath.arrayAddElement(strings, string);
    }

    static public String[] stringListAddList(String[] strings, String[] newStrings)
    {
        return (String[])StsMath.arrayAddArray(strings, newStrings);
    }
    static public boolean stringListHasString(String[] strings, String string)
    {
        if(strings == null || string == null) return false;
        for(int n = 0; n < strings.length; n++)
            if(strings[n].equals(string)) return true;
        return false;
    }

    static public String trimSuffix(String string, String[] suffixes)
    {
        for(int n = 0; n < suffixes.length; n++)
        {
            int index = string.toLowerCase().indexOf(suffixes[n].toLowerCase());
            if(index > 0) return string.substring(0, index);
        }
        return null;
    }


    static public String trimSuffix(String string)
    {
        int dotIndex = string.lastIndexOf(".");
        if(dotIndex == -1) return string;
        return string.substring(0, dotIndex);
    }

    static public String trimPrefix(String string, String prefix)
    {
        if(!string.startsWith(prefix)) return string;
        int prefixLength = prefix.length();
        int length = string.length();
        return string.substring(prefixLength, length);
    }

    static public boolean stringContainsString(String string, String subString)
    {
		if(string == null || subString == null) return false;
        return string.indexOf(subString) >= 0;
    }

    static public String padClipString(String string, int num)
    {
        if(string.length() > num)
            string = string.substring(0,num);
        else
        {
            for (int i = string.length(); i < num; i++)
                string = string + " ";
        }
        return string;
    }

	static public void addFixedLengthStringX(StringBuffer stringBuffer, String string, int num)
	 {
		 char[] chars = new char[num];
		 char[] stringChars = string.toCharArray();
		 int length = Math.min(num, stringChars.length);
		 System.arraycopy(stringChars, 0, chars, 0, length);
		 stringBuffer.append(chars);
    }

	static public void addFixedLengthString(StringBuffer stringBuffer, String string, int num)
	 {
		 if(string.length() > num)
		 {
			 string = string.substring(0, num);
			 stringBuffer.append(string);
		 }
		 else
		 {
			 StringBuffer newStringBuffer = new StringBuffer(num);
			 newStringBuffer.append(string);
			 for(int n = string.length(); n < num; n++)
				 newStringBuffer.append(' ');
			 stringBuffer.append(newStringBuffer);
		 }
    }
	static public String matchStringInList(String[] strings, String string, String nullString)
	{
		if(string == null) return nullString;
		for(int n =0; n < strings.length; n++)
			if(strings[n].equals(string)) return strings[n];
		return nullString;
	}
	static public final boolean stringsEqual(String s1, String s2)
	{
		return s1 == s2 || s1.equals(s2);
	}

    static public boolean isNumeric(String token)
    {
        String numericValues = new String("0123456789.,E");
        char[] oneChar = token.toCharArray();
        for(int i=0; i<token.length(); i++)
        {
            if(numericValues.indexOf(oneChar[i]) < 0)
                return false;
        }
        return true;
    }
    static public String[] getTokens(String line)
    {
        if(line == null) return null;
        line.trim();
        int nTokens = 0;
        StringTokenizer stringTokenizer =  new StringTokenizer(line);
        nTokens = stringTokenizer.countTokens();
        if(nTokens == 0) return null;
        String[] tokens = new String[nTokens];
        int n = 0;
        while(stringTokenizer.hasMoreTokens())
            tokens[n++] = stringTokenizer.nextToken();
        return tokens;
    }

    static public String[] getTokens(String line, String delimiter)
    {
        if(line == null) return null;
        line.trim();
        int nTokens = 0;
        StringTokenizer stringTokenizer =  new StringTokenizer(line, delimiter);
        nTokens = stringTokenizer.countTokens();
        if(nTokens == 0) return null;

        String[] tokens = new String[nTokens];
        int n = 0;
        while(stringTokenizer.hasMoreTokens())
            tokens[n++] = stringTokenizer.nextToken();
        return tokens;
    }

    static public String[] getTokens(String line, StringTokenizer stringTokenizer)
    {
        if(line == null) return null;
        line.trim();
        int nTokens = 0;
        nTokens = stringTokenizer.countTokens();
        if(nTokens == 0) return null;
        String[] tokens = new String[nTokens];
        int n = 0;
        while(stringTokenizer.hasMoreTokens())
            tokens[n++] = stringTokenizer.nextToken();
        return tokens;
    }
    static public String deWhiteSpaceString(String s)
    {
        s = detabString(s);
        return deSpaceString(s);
    }

    static public String cleanLine(String s, String delimiters)
    {
        if(s == null) return null;
        String res = "";
        int len = s.length();
        int pos = 0;
        boolean foundFirst = false;
        for (int i=0; i < len; i++)
        {
            char c = s.charAt(i);
            if(delimiters.contains(s.substring(i,i+1)))
            {
                if(foundFirst)
                    res  += " " + c;
                else
                {
                    foundFirst = true;
                    res += c;
                }
            }
            else
            {
                res += c;
                foundFirst = false;
            }
        }
        if(foundFirst)
            res += " ";
        return res;
    }

    static public String detabString(String s)
    {
        if(s == null) return null;
        if (s.indexOf('\t') == -1)
            return s;
        String res = "";
        int len = s.length();
        int pos = 0;
        int i = 0;
        for (; i < len && s.charAt(i) == '\t'; i++)
        {
            res += "        ";
            pos += 8;
        }
        for (; i < len; i++)
        {
            char c = s.charAt(i);
            if (c == '\t') {
                do {
                    res += " ";
                    pos++;
                } while (pos % 8 != 0);
            }
            else {
                res += c;
                pos++;
            }
        }
        return res;
     }

    static public String deQuoteString(String s)
    {
        if(s == null) return null;
        if (s.indexOf('"') == -1)
            return s;
        String res = "";
        int len = s.length();
        int pos = 0;
        int i = 0;
        int start = 0, end = len-1;
        if(s.charAt(0) == '"')
            start = 1;
        if(s.charAt(end) == '"')
            end = end - 1;
        for(int ii=start; ii<=end; ii++)
            res += s.charAt(ii);
        return res;
     }

    static public String deSpaceString(String s)
    {
        if(s == null) return null;
        if (s.indexOf(' ') == -1)
            return s;
        String res = "";
        int len = s.length();
        int pos = 0;
        int i = 0;
        int start = 0, end = len-1;
        if(s.charAt(0) == ' ')
            start = 1;
        if(s.charAt(end) == ' ')
            end = end - 1;
        for(int ii=start; ii<=end; ii++)
            res += s.charAt(ii);
        return res;
     }

    static public int compareTokenStrings(String string1, String string2, String delimiter)
    {
        String[] tokens1 = getTokens(string1, delimiter);
        String[] tokens2 = getTokens(string2, delimiter);
        int nTokens = Math.min(tokens1.length, tokens2.length);
        for(int n = 0; n < nTokens; n++)
        {
            int compare = compareStrings(tokens1[n], tokens2[n]);
            if(compare != 0) return compare;
        }
        return tokens1.length - tokens2.length;
    }

    static public int compareStrings(String s1, String s2)
    {
        double d1 = getStringDouble(s1);
        if(d1 == StsParameters.largeDouble)
            return s1.compareTo(s2);
        double d2 = getStringDouble(s2);
        if(d2 == StsParameters.largeDouble)
            return s1.compareTo(s2);
        return (int)(d1 - d2);
    }

    static public double getStringDouble(String string)
    {
        try { return Double.parseDouble(string); }
        catch(Exception e) { return StsParameters.largeDouble; }
    }

    static public void main(String[] args)
	{
        int compare;
        String s1 = "3.10.0";
        String s2 = "3.10";
        testCompare(s1, s2, ".");
        s2 = "3.9.3";
        testCompare(s1, s2, ".");
        s2 = "3.10-tom";
        testCompare(s1, s2, ".");
        /*
        String string = "file:abc";
        String newString = StsStringUtils.trimPrefix(string, "file:");
        System.out.println("new: " + newString);
        */
    /*
        StringBuffer stringBuffer = new StringBuffer(3200);
		addFixedLengthString(stringBuffer, "C1", 80);
		addFixedLengthString(stringBuffer, "123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789X", 80);
		String string = stringBuffer.toString();
		System.out.println(string);
		System.out.println("string length: " + string.length());
    */
	}

    static private void testCompare(String s1, String s2, String delimiter)
    {
       int compare = StsStringUtils.compareTokenStrings(s1, s2, delimiter);
       if(compare < 0)
            System.out.println(s1 + " is older than " + s2);
       else if(compare > 0)
            System.out.println(s1 + " is newer than " + s2);
       else
           System.out.println(s1 + " is same as " +  s2);
    }
}
