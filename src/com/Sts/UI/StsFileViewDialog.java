package com.Sts.UI;

import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.font.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsFileViewDialog extends JDialog {
    JPanel panel1 = new JPanel();
    JLabel fileViewLbl = new JLabel();
    JTextArea fileTextArea = new JTextArea();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    JScrollPane jScrollPane1 = new JScrollPane();
    int height = 500;
    int width = 500;
    Font font = new java.awt.Font("Serif", 1, 12);
    LineMetrics fm = null;

    public StsFileViewDialog(Frame frame, String title, boolean modal)
    {
        this(frame, title, modal, 500, 500);
    }
    public StsFileViewDialog(Frame frame, String title, boolean modal, int width, int height)
    {
        super(frame, title, modal);
        try
        {
            this.width = width;
            this.height = height;
            jbInit();
            pack();
			StsToolkit.centerComponentOnScreen(this);
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public StsFileViewDialog()
    {
        this(null, "", false);
    }

    private void jbInit() throws Exception
    {
        panel1.setLayout(gridBagLayout1);
        fileViewLbl.setFont(new java.awt.Font("Serif", 1, 12));
        fileViewLbl.setHorizontalAlignment(SwingConstants.CENTER);
        fileViewLbl.setText("File View");
        fileTextArea.setBackground(Color.lightGray);
        fileTextArea.setEditable(false);
        fileTextArea.setText("");
        fileViewLbl.setFont(font);
        fileTextArea.setLineWrap(false);
        fileTextArea.setFont(new Font("Monospaced", 1, 11));
    fileTextArea.setBorder(BorderFactory.createEtchedBorder());
    fileTextArea.setDoubleBuffered(true);
    fileTextArea.setCaretColor(Color.red);
    fileTextArea.setCaretPosition(0);
        jScrollPane1.setPreferredSize(new Dimension(width,height));

        jScrollPane1.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        jScrollPane1.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        jScrollPane1.getViewport().add(fileTextArea);
        getContentPane().add(panel1, BorderLayout.CENTER);
        panel1.add(fileViewLbl,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        panel1.add(jScrollPane1,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(8, 6, 9, 7), 0, 0));
    }

    public void setViewTitle(String Title)
    {
        fileViewLbl.setText(Title);
    }

    public void setFileText(String headerText)
    {
        fileTextArea.setText(headerText);
    }

    public void appendLine(String line)
    {
        fileTextArea.append(line + '\n');
    }

    public void setEditable(boolean edit)
    {
        if(edit)
        {
            fileTextArea.setEditable(true);
            fileTextArea.setBackground(SystemColor.white);
        }
        else
        {
            fileTextArea.setBackground(Color.lightGray);
            fileTextArea.setEditable(false);
        }
    }
    public String getText()
    {
        return fileTextArea.getText();
    }
}
