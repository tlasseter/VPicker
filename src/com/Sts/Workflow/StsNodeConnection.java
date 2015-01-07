package com.Sts.Workflow;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

/** provides the connection between a bundle of input nodes and an output node
 *
 */

public class StsNodeConnection
{
    StsNodeBundle inNodes;
    StsWorkflowTreeNode outNode;

    public StsNodeConnection(StsWorkflowTreeNode inNode, StsWorkflowTreeNode outNode)
    {
        inNodes = new StsNodeBundle(inNode);
        this.outNode = outNode;
        inNodes.addOutputNode(outNode);
        outNode.addInputNodeBundle(inNodes);
    }

    public StsNodeConnection(StsNodeBundle inNodes, StsWorkflowTreeNode outNode)
    {
        this.inNodes = inNodes;
        this.outNode = outNode;
        inNodes.addOutputNode(outNode);
        outNode.addInputNodeBundle(inNodes);
    }

    protected boolean adjustWorkflowPanelState(String actionClassName, int actionStatus)
    {
        boolean changed = false;
        if(inNodes.adjustedState(actionClassName, actionStatus))
        {
            if(outNode.updateActionStatus()) changed = true;
        }
        else if(outNode.adjustedState(actionClassName, actionStatus))
        {
            ;
// This is not needed anymore. It would set all pre-requisites to ended
//            if(outNode.actionStatus > StsWorkflowTreeNode.STARTED)
//            {
//                inNodes.setActionStatus(StsWorkflowTreeNode.ENDED);
//                changed = true;
//            }

        }
        return changed;
    }
}