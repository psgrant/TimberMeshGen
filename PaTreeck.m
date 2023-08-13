classdef PaTreeck
    properties
        % Data that is stored at each node
        Data
        % List of each nodes parent
        Parent
        % To see if the branch is complete 1 for complete
        BranchComplete
        % What level each node is - useful for plotting
        Level
        % X Position of node in the tree
        XPos
    end
    
    % Start class methods
    methods
        % Create tree and initalise with the root data
        function [T, RootID] = PaTreeck(Data)
            T.Parent = 0;
            T.Data = cell(1);
            T.Data{1} = Data;
            RootID = 1;
            T.BranchComplete = 0;
        end
        
        % Add a new node (leaf) to the tree
        function [T, NodeID] = addNode(T,Data, Parent)
            
            % Update the Parent ID list
            T.Parent = [T.Parent; Parent];
            % Update with an entry in the branch end completed vector
            T.BranchComplete = [T.BranchComplete; 0];
            % Input the data into the tree
            T.Data{end + 1,1} = Data;
            % Return the number assigned to that node
            NodeID = numel(T.Data);         
        end
        
        % return the number of children a node has
        function [Parent] = getParent(T,NodeID)
            % return the parent of the node
            Parent = T.Parent(NodeID);
        end
        
        % return the NodeIds of the children of a branch. Will return an
        % empty column vector if the node has no children
        function [ChildrenIDs] = getChildren(T, NodeID)
            % Find the children Ids of the node
            ChildrenIDs = find(T.Parent == NodeID);
        end
        
        % get function returning the data stored at the given node
        function [Content] = get(T,NodeID)
            
            Content = T.Data{NodeID};
        end
        
        % set function to set the content at the given node
        function [T] = set(T, Content, NodeID)
            % Set the content type
            T.Data{NodeID} = Content;
        end
        
        % Returns the list of siblings (not inclduing the starting node)
        % will return an empty colum vector if the node has no siblings
        function [SiblingList] = getSiblings(T,NodeID)
            % Find the siblings (children of the aprent node)
            SiblingList = T.getChildren(T.getParent(NodeID));
            % Remove the NodeId of the input node
            SiblingList(SiblingList == NodeID) = [];
        end
        
        % Updates the branch complete property of the tree. It also checks
        % to see if the sibling node is at the end and if so will mark the
        % parent as complete, also returning the next node that is in line
        % in the tree
        function [T, NextNode] = endOfBranch(T,NodeID)
            T.BranchComplete(NodeID) = 1;
            FoundNewNode = 0;
            % Get the sibling node
            SiblingList = T.getSiblings(NodeID);
            
            % If the sibling node has not been completed (or checked for
            % completeness)
            if T.BranchComplete(SiblingList) == 0
                NextNode = SiblingList;
                
            else % The sibling branch has been completed
                
                % Update the parent node
                T.BranchComplete(T.Parent(NodeID)) = 1;
                % The Next node is the sibling of the parent node
                NextNode = T.getSiblings(T.getParent(NodeID));
                
                % We need to loop through the tree until we find a node
                % that has not been completed
                while FoundNewNode == 0
                    
                    % If both siblings have a completed branch we move up the
                    % tree until we find a valid next node or we complete the
                    % tree
                    if T.BranchComplete(NextNode) == 1
                        % Update the parent node to have a complete branch
                        T.BranchComplete(T.Parent(NextNode)) = 1;
                        
                        NextNode = T.getSiblings(T.getParent(NextNode));
                    else
                        % The new node is not complete (or checked for
                        % completeness) break out of the loop and we have
                        % our new node
                        FoundNewNode = 1;
                        
                    end
                end
            end
        end
    end
    
    
    
end

