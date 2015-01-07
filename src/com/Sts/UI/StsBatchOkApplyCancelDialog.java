package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;

import java.awt.*;

/** This is a dialog centered on a parent with a displayPanel displayed and an okCancelObject which it communicates with
 *  when ok/apply/cancel buttons are pushed.
 */

public class StsBatchOkApplyCancelDialog extends StsOkApplyCancelDialog
{
    private StsButton submitButton = new StsButton("Batch", "Save settings and start batch process of data", this, "submit");
    private boolean batchProcess = false;
    private StsPreStackLineSet volume;
    private byte batchType = -1;

    public StsBatchOkApplyCancelDialog(Frame frame, StsPreStackLineSet volume, byte batchType, StsDialogFace[] okCancelObjects, String title, boolean modal)
    {
        super(frame, title, modal);
        if(frame != null) setLocation(frame.getLocation());
        this.volume = volume;
        this.batchType = batchType;
        this.okCancelObjects = okCancelObjects;
        layoutPanels();
		super.addWindowCloseOperation();
		display();
    }

    private void layoutPanels()
    {
        panel.gbc.fill = GridBagConstraints.HORIZONTAL;
        for(int n = 0; n < okCancelObjects.length; n++)
        {
            okCancelObjects[n] = okCancelObjects[n].getEditableCopy();
            panel.add(okCancelObjects[n].getPanel());
        }
        Insets insets = new Insets(4, 4, 4, 4);
        okButton.setMargin(insets);
        cancelButton.setMargin(insets);
        applyButton.setMargin(insets);
        submitButton.setMargin(insets);

        buttonBox.addToRow(okButton);
        buttonBox.addToRow(applyButton);
        buttonBox.addToRow(cancelButton);
        buttonBox.addEndRow(submitButton);
        panel.add(buttonBox);

        getContentPane().add(panel);
	}

    public void submit()
    {
        if(StsYesNoDialog.questionValue(this,"This will take some time. Do you wish to continue?"))
        {
            // Save all the current settings
            for(int n = 0; n < okCancelObjects.length; n++)
                okCancelObjects[n].dialogSelectionType(StsDialogFace.OK);

            // Dispose of the dialog
            dispose();
            volume.batchProcess(batchType);
        }
	}

    public boolean batchProcess() { return batchProcess; }
}
