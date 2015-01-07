package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version 1.0
 */

public class StsOverwriteDialog extends JDialog
{
	public static final int SKIP_FILE = 0;
	public static final int USE_FILE = 1;

	private StsJPanel panel = StsJPanel.addInsets();
	private JLabel lineOne;
    private JLabel lineTwo;
    private StsButtonListFieldBean useKeepRadioButtons;
    private StsBooleanFieldBean useOptionCheckbox;
    private StsButton okButton = new StsButton("OK", "Accept changes and dismiss dialog.", this, "ok");

    private String filename;
    String removeKeepString =  removeKeepStrings[0];
    private boolean useOption = false;

    static private String[] removeKeepStrings = new String[] { "Don't use this file (remove from list)", "Keep file" };
    public StsOverwriteDialog(Frame owner, String title, boolean modal, String filename)
	{
		super(owner, title, modal);
		try
		{
			this.filename = filename;
			buildPanel();
            pack();
        }
		catch (Exception exception)
		{
			exception.printStackTrace();
		}
	}

	public StsOverwriteDialog(String filename)
	{
		this(null, "Already Processed Files Dialog", true, filename);
	}

	public void setUseOption(boolean use)
	{
		this.useOption = use;
	}

    public boolean getUseOption() { return useOption; }

    public void setRemoveKeepString(String string) { removeKeepString = string; }
    public String getRemoveKeepString() { return removeKeepString; }
    public int getResponse()
    {
        for(int n = 0; n < 2; n++)
            if(removeKeepString == removeKeepStrings[n])
                return n;
        return 0;
    }

    public boolean getDontAskAgain() { return useOption; }
    
    private void buildPanel() throws Exception
	{
		this.setSize(350, 200);
		this.setResizable(false);
		this.setLocationRelativeTo(null);
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        panel.gbc.anchor = GridBagConstraints.WEST;
        lineOne = new JLabel(filename + " has already been processed.");
        lineOne.setHorizontalTextPosition(SwingConstants.LEFT);
        lineTwo = new JLabel("What would you like to do?");
        lineTwo.setHorizontalTextPosition(SwingConstants.LEFT);
        panel.add(lineOne);
        panel.add(lineTwo);

        panel.gbc.anchor = GridBagConstraints.CENTER;
        panel.gbc.insets = new Insets(0, 20, 0, 0);
        useKeepRadioButtons = new StsButtonListFieldBean(this, "removeKeepString", null, removeKeepStrings, false);
        panel.add(useKeepRadioButtons);

        panel.gbc.anchor = GridBagConstraints.WEST;
        panel.gbc.insets = new Insets(0, 0, 0, 0);
        useOptionCheckbox = new StsBooleanFieldBean(this, "useOption", "Always use this option (for remainder of this session", false);
        panel.add(useOptionCheckbox);

        Insets insets = new Insets(4, 4, 4, 4);
		okButton.setMargin(insets);
        panel.gbc.anchor = GridBagConstraints.CENTER;
        panel.add(okButton);

        add(panel);
    }

	public void ok()
	{
		dispose();
	}

    static public void main(String[] args)
    {
        try
        {
            UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
            //UIManager.setLookAndFeel(new com.sun.java.swing.plaf.motif.MotifLookAndFeel());
            //UIManager.setLookAndFeel(new javax.swing.plaf.metal.MetalLookAndFeel());
        }
        catch (Exception e)
        {
        }

        Frame frame = new Frame("Test Frame");
        frame.setSize(new Dimension(600, 400));
        StsToolkit.centerComponentOnScreen(frame);
        StsOverwriteDialog dialog = new StsOverwriteDialog(frame, "Test", false, "filename");
        dialog.setVisible(true);
    }
}
