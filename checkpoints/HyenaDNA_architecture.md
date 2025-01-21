## Embedding layer
- [GPT2Embeddings class](/for_testing.ipynb#C17:L6)

    The embedding layer has implemented intermediate projections to a **word_embe_proj_dim dimension** space (in case that parameter is not None) and then to the final **embed_dim dimension** (d_model). Presumably this allows to handle pre-trained embeddings that have a different dimension than the one handled by HyenaDNA.
    HyenaDNA implement positional encoding within the Hyena operator, but the embedding layer allows it if the **max_position_embeddings** parameter > 0. In that case it uses positional embedding for spatial encoding.

## Backbone
The backbone of the HyenaDNA addmite two ways of normalization: pre-normalization or post-normalization. Pre- and post- refers to the position of de normalization layer with respect to mixer layer. The code has the flexibility to use a Hyena operator or multiead attention (MHA) instead as the mixer layer. As pre-normalization seems to have a better performance, here we'll describe that architecture.

The backbone of the model consists of an embbeding layer fallowed by a stack of blocks with the structure we describe below. 
- [Backbone](/for_testing.ipynb#C17:L252)

    ::: mermaid
    flowchart LR
        input[input]:::noBorder
        hidden_state_out[hidden_state]:::noBorder
        plus(("&plus;"))
        Embedding[Embedding]
        Layer1[Block]
        Layer2[Block]
        Layer_dots["..."]:::noBorder
        Layer_last[Block]:::BlockSty
        Norm[Norm]

        input --> Embedding
        Embedding -- hidden_state --> Layer1
        Layer1 -- residual --> Layer2
        Layer1 -- hidden_state --> Layer2
        Layer2 --> Layer_dots
        Layer2 --> Layer_dots
        Layer_dots --> Layer_last
        Layer_dots --> Layer_last
        Layer_last --> plus
        Layer_last --> plus
        plus --> Norm
        Norm --> hidden_state_out

    classDef noBorder stroke-width:0px,fill:none;
    :::

- [Block](/for_testing.ipynb#C17:L62)

    (**prenorm=True**, **mixer_subset = None**, and ignoring dropout layers)

    ::: mermaid
    flowchart LR
        input[input]:::noBorder
        residual_out["residual"]:::noBorder
        hidden_state_out["hidden_state"]:::noBorder
        plus(("&plus;"))
        Mixer[Hyena]
        Norm1[Norm]
        Norm2[Norm]

        subgraph outputs [ ]
            direction TB
            residual_out
            hidden_state_out
        end 
        
        input -- residual --> plus
        input -- hidden_state --> Norm1
        Norm1 --> Mixer
        Mixer --> plus
        plus --> residual_out
        plus --> Norm2
        Norm2 --> MLP
        MLP --> hidden_state_out

    classDef noBorder stroke-width:0px,fill:none;
    style outputs fill:none,stroke:none;
    :::



- [Hyena operator](/for_testing.ipynb#C17:L141-L215)

    Here **d_model** is the emmbedding dimension. The inner vectors **x0**, **x1**, **v** had the same **d_model** dimension.
    (**inner_width = d_model * (order + 1)**, **order = 2**)

    ::: mermaid
    flowchart LR
        input[input]:::noBorder
        output[output]:::noBorder
        product1(("&times;"))
        product2(("&times;"))
        Conv(Conv_1D):::vertical
        Filter["Hyena
                Filter"]
        in_proj@{ shape: sl-rect, label: "Dense
                d_model x inner_width" }
        out_proj["Dense
                d_model x d_model"]

        subgraph split["split
                        inner vector"]
            direction TB
            x0:::noBorder
            x1:::noBorder
            v:::noBorder
        end 
        
        input --> in_proj
        in_proj --> Conv
        Conv --> split
        x0 --> product2
        x1 --> product1
        v --> product1
        product1 --> Filter
        Filter --> product2
        product2 --> out_proj
        out_proj --> output

    classDef noBorder stroke-width:0px,fill:none;
    classDef vertical height:130px, label-rotate: 90, display:flex, align-items:center, justify-content:center;
    style split fill:none, display:flex, vertical-align: middle;
    :::

- [Hyena Filter](/for_testing.ipynb#C17:L63-L139)

    ::: mermaid
    flowchart LR
        input[input]:::noBorder
        output[output]:::noBorder
        t:::noBorder
        
        Pos_Emb(Positional
                Embedding)
        Filter["implicit filter
                (FFN)"]
        Modulation["z exp(-a t)
                modulation"]
        fftconv[FFT
                convolution]

        t --> Modulation
        input -->  fftconv
        Pos_Emb -- z --> Filter
        Filter -- z --> Modulation
        Modulation -- k --> fftconv
        fftconv --> output
        

    classDef noBorder stroke-width:0px,fill:none;
    :::

