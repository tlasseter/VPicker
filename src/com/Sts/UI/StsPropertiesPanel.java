package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public class StsPropertiesPanel extends StsJPanel
{
	StsButton expandButton;
	StsButton dismissButton = new StsButton("Close", StsToolkit.getCloseIcon(), "Close panel.", this, "dismissPanel");
	StsGroupBox panelWithDismiss;
    int numberOfBeans = 0;
	
	public StsPropertiesPanel(StsFieldBean[] beans, String title)
	{
		this(beans, title, false);
	}
	public StsPropertiesPanel(StsFieldBean[] beans, String title, boolean expandIt)
	{
		dismissButton.setPreferredSize(new Dimension(20, 20));
		expandButton = new StsButton(title, "Select to expand panel.", this, "expandPanel");
		expandButton.setMargin(new Insets(2, 2, 2, 2));
		panelWithDismiss = new StsGroupBox(title);
		panelWithDismiss.gbc.anchor = GridBagConstraints.WEST;
		panelWithDismiss.add(dismissButton);
		StsJPanel panel = StsJPanel.addInsets();
		panel.gbc.anchor = GridBagConstraints.WEST;
        panel.gbc.fill = GridBagConstraints.HORIZONTAL;
		panel.add(beans);
        numberOfBeans = beans.length;
        gbc.fill = gbc.HORIZONTAL;
		panelWithDismiss.add(panel);
		add(expandButton);
		if(expandIt) expandPanel();
	}

	public void expandPanel()
	{
        collapseAll();
		removeAll();
        remove(expandButton);

        // Required so the properties panel doesn't expand beyond the bottom of the
        // screen or beyond the vertical dimensions of the screen
        JScrollPane scrollPane = new JScrollPane();
		scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
		scrollPane.getViewport().add(panelWithDismiss);
        if(numberOfBeans > 12)
            scrollPane.setPreferredSize(new Dimension(200, 400));
		add(scrollPane);

		rebuild();
	}

	private void rebuild()
	{
		Component parent = this;
		while(!(parent instanceof JDialog))
		{
			parent = parent.getParent();
			if(parent == null) return;
		}
		JDialog dialog = (JDialog)parent;
		dialog.pack();
	}

    private void collapseAll()
    {
        Component parent = this;
        while(!(parent instanceof StsOkApplyCancelDialog))
        {
            parent  = parent.getParent();
            if(parent == null) return;
        }
        ((StsOkApplyCancelDialog)parent).collapseAll();
    }

	public void dismissPanel()
	{
		removeAll();
		add(expandButton);
		rebuild();
	}

	public static void main(String[] args)
	{
		try { UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel()); }
	    catch(Exception e) { }
		StsModel model = new StsModel();
		model.createSpectrums();
		StsProject project = new StsProject();
		model.setProject(project);
		StsPreStackLineSet3dClass seismicClass = new StsPreStackLineSet3dClass();
		seismicClass.projectInitialize(model);
		StsSemblanceDisplayProperties semblanceDisplayProperties = new StsSemblanceDisplayProperties(StsPreStackLineSet3d.DISPLAY_MODE, null);
		StsWiggleDisplayProperties wiggleProperties = new StsWiggleDisplayProperties(seismicClass, null);
		StsFilterProperties filterProperties = new StsFilterProperties(null);
		StsOkCancelDialog okCancelDialog = new StsOkCancelDialog(null, new StsDialogFace[] {semblanceDisplayProperties, wiggleProperties, filterProperties}, "Semblance Properties Test", false);
	}
}
