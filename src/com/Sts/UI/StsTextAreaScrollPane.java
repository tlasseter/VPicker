package com.Sts.UI;

import com.Sts.Utilities.*;

import javax.swing.*;
import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jan 15, 2010
 * Time: 7:46:16 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsTextAreaScrollPane extends JScrollPane
{
    public StsTextArea textArea;

    public StsTextAreaScrollPane()
    {
        textArea = new StsTextArea();
        setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
		getViewport().add(textArea);
    }

    public void append(String string)
    {
        textArea.append(string);
    }

    public void read(FileReader fileReader, Object object)
    {
        try
        {
            textArea.read(fileReader, object);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "read", e);
        }
    }
}
