package com.Sts.WorkflowPlugIn.Html;

import com.Sts.Utilities.*;

import java.net.*;

/**
 * Title:        Workflow development
 * Description:
 * Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */

public class StsHtml
{
    static StsHtml htmlInstance = new StsHtml();
    static Class htmlClass = htmlInstance.getClass();

    public StsHtml()
    {
    }
    static public URL getHtml(String name)
    {
        return getHtml(htmlClass, name);
    }

    static public URL getHtml(Class c, String filename)
    {
        try {
            java.net.URL url = c.getResource(filename);
            if (url == null)
                throw new StsException();
            return url;
        }
        catch (Exception e) {
            StsException.systemError("StsHtml.getHtml() failed.\n" +
                                     "Html file not found: " + filename);
            return null;
        }
    }
}