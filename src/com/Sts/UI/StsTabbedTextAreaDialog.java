package com.Sts.UI;

import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.io.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsTabbedTextAreaDialog extends JDialog
{
	private static final int DEFAULT_X_DIMENSION = 500;
	private static final int DEFAULT_Y_DIMENSION = 500;

	private JPanel panel1 = new JPanel();
	private JTabbedPane tabbedPane = new JTabbedPane();
	private BorderLayout borderLayout1 = new BorderLayout();

	public StsTabbedTextAreaDialog(Frame frame, String title, boolean modal)
	{
		this(frame, title, modal, DEFAULT_X_DIMENSION, DEFAULT_Y_DIMENSION);
	}

	public StsTabbedTextAreaDialog(Frame frame, String title, boolean modal, int width, int height)
	{
		super(frame, title, modal);
		try
		{
			this.setSize(new Dimension(width, height));
			jbInit();
			StsToolkit.centerComponentOnScreen(this);
            setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

	public StsTabbedTextAreaDialog()
	{
		this(null, "", false);
	}

    public StsTabbedTextAreaDialog(Frame frame, String title, boolean modal, int width, int height, File[] files)
	{
        this(frame, title, modal, width, height);
        for(int n = 0; n < files.length; n++)
            addFile(files[n]);
    }

	private void jbInit() throws Exception
	{
		panel1.setLayout(borderLayout1);
		getContentPane().add(panel1, BorderLayout.CENTER);
		panel1.add(tabbedPane, java.awt.BorderLayout.CENTER);
	}

	public void setTabTitle(int tabIndex, String title)
	{
		tabbedPane.setTitleAt(tabIndex, title);
	}

	public void setText(int tabIndex, String headerText)
	{
		getTextArea(tabIndex).setText(headerText);
	}

	public void appendLine(int tabIndex, String line)
	{
		getTextArea(tabIndex).append(line + '\n');
	}

	public void appendLine(JTextArea textArea, String line)
	{
		textArea.append(line + '\n');
	}

    public void addFile(File file)
    {
        if(file == null) return;
        StsTextAreaScrollPane textArea = addTab(file.getName());
        try
        {
            FileReader fileReader = new FileReader(file);
            textArea.read(fileReader, null);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "readFile", e);
        }
    }
	public void setEditable(int tabIndex, boolean edit)
	{
		JTextArea textArea = getTextArea(tabIndex);
		if (edit)
		{
			textArea.setEditable(true);
			textArea.setBackground(SystemColor.white);
		}
		else
		{
			textArea.setBackground(Color.lightGray);
			textArea.setEditable(false);
		}
	}

	public String getText(int tabIndex)
	{
		return getTextArea(tabIndex).getText();
	}

	public StsTextAreaScrollPane addTab(String title)
	{
		StsTextAreaScrollPane textAreaScrollPane = new StsTextAreaScrollPane();
		tabbedPane.addTab(title, textAreaScrollPane);
        return textAreaScrollPane;
    }

	public int getTabIndex(String title)
	{
		return tabbedPane.indexOfTab(title);
	}

	public void setTabActive(int tabIndex)
	{
		tabbedPane.setSelectedIndex(tabIndex);
	}

	public JTextArea getTextArea(int tabIndex)
	{
		return (JTextArea)((JScrollPane)tabbedPane.getComponentAt(tabIndex)).getViewport().getView();
	}

    static public void main(String[] args)
    {
        StsFileChooser fileChooser = StsFileChooser.createFileChooserPostfix(null, "pick text files", ".", ".txt");
        fileChooser.show();
        File[] files = fileChooser.getFiles();
        int nFiles = files.length;
        StsTabbedTextAreaDialog tabbedDialog = new StsTabbedTextAreaDialog(null, "Test Tabbed Text Area Dialog", true, 500, 500, files);
        tabbedDialog.setVisible(true);
    }
}
