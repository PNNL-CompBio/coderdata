import torch
import torch.nn as nn
import torch.nn.functional as F
import pickle


class MultiheadAttention(nn.Module):

    def __init__(self, input_dim, embed_dim, num_heads):
        super().__init__()
        assert embed_dim % num_heads == 0, "Embedding dimension must be 0 modulo number of heads."

        self.embed_dim = embed_dim # 512
        self.num_heads = num_heads # 8
        self.head_dim = embed_dim // num_heads

        # Stack all weight matrices 1...h together for efficiency
        # Note that in many implementations you see "bias=False" which is optional
        self.qkv_proj = nn.Linear(input_dim, 3*embed_dim)
        self.o_proj = nn.Linear(embed_dim, embed_dim)

        self._reset_parameters()

    def _reset_parameters(self):
        # Original Transformer initialization, see PyTorch documentation
        nn.init.xavier_uniform_(self.qkv_proj.weight)
        self.qkv_proj.bias.data.fill_(0)
        nn.init.xavier_uniform_(self.o_proj.weight)
        self.o_proj.bias.data.fill_(0)

    def forward(self, x, mask=None, attn_bias=None, return_attention=False):
        batch_size, seq_length, _ = x.size()
        # mask = x.eq(0)
        # if mask is not None:
        #     mask = expand_mask(mask)
        x = x * (1 - mask.unsqueeze(-1).type_as(x))
        qkv = self.qkv_proj(x)

        # Separate Q, K, V from linear output
        qkv = qkv.reshape(batch_size, seq_length, self.num_heads, 3*self.head_dim)
        qkv = qkv.permute(0, 2, 1, 3) # [Batch, Head, SeqLen, Dims]
        q, k, v = qkv.chunk(3, dim=-1)

        # Determine value outputs
        # values, attention = scaled_dot_product(q, k, v, mask=mask)
        d_k = q.size()[-1]
        attn_logits = torch.matmul(q, k.transpose(-2, -1))
        attn_logits = attn_logits * ( d_k**(-0.5) )
        attn_logits = attn_logits.masked_fill(mask.unsqueeze(1).unsqueeze(2).to(torch.bool), float('-inf') )
        attn_logits = attn_logits.view(batch_size*self.num_heads, seq_length, seq_length)
        
        
        # attention_mask = torch.randn(4,4,4)
        
        if attn_bias:
            attn_bias = attn_bias.unsqueeze(1)
            attn_bias = attn_bias.repeat(1, self.num_heads, seq_length, seq_length)
            attn_bias = attn_bias.view(batch_size*self.num_heads, seq_length, seq_length)
            attn_logits = attn_logits + attn_bias

        attn_logits = attn_logits.view(batch_size, self.num_heads, seq_length, seq_length)
        attention = F.softmax(attn_logits, dim=-1)
        
        values = torch.matmul(attention, v)
        values = values.permute(0, 2, 1, 3) # [Batch, SeqLen, Head, Dims]
        values = values.reshape(batch_size, seq_length, self.embed_dim)

        
        # values = values.permute(0, 2, 1, 3) # [Batch, SeqLen, Head, Dims]
        # values = values.reshape(batch_size, seq_length, self.embed_dim)
        o = self.o_proj(values)

        if return_attention:
            return o, attention
        else:
            return o
        
        
class EncoderBlock(nn.Module):

    def __init__(self, input_dim, num_heads, dim_feedforward, dropout=0.0):
        """
        Inputs:
            input_dim - Dimensionality of the input
            num_heads - Number of heads to use in the attention block
            dim_feedforward - Dimensionality of the hidden layer in the MLP
            dropout - Dropout probability to use in the dropout layers
        """
        super().__init__()

        # Attention layer
        self.self_attn = MultiheadAttention(input_dim, input_dim, num_heads)

        # Two-layer MLP
        self.linear_net = nn.Sequential(
            nn.Linear(input_dim, dim_feedforward),
            nn.Dropout(dropout),
            nn.ReLU(inplace=True),
            nn.Linear(dim_feedforward, input_dim)
        )

        # Layers to apply in between the main layers
        self.norm1 = nn.LayerNorm(input_dim)
        self.norm2 = nn.LayerNorm(input_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, mask=None, attn_bias=None):
        # Attention part
        attn_out = self.self_attn(x, mask=mask, attn_bias=attn_bias)
        x = x + self.dropout(attn_out)
        x = self.norm1(x)

        # MLP part
        linear_out = self.linear_net(x)
        x = x + self.dropout(linear_out)
        x = self.norm2(x)

        return x
    
    
class TransformerEncoder(nn.Module):

    def __init__(self, num_layers, **block_args):
        super().__init__()
        self.layers = nn.ModuleList([EncoderBlock(**block_args) for _ in range(num_layers)])

    def forward(self, x, mask=None, attn_bias=None):
        for l in self.layers:
            x = l(x, mask=mask, attn_bias=attn_bias)
        return x

    def get_attention_maps(self, x, mask=None, attn_bias=None):
        attention_maps = []
        for l in self.layers:
            _, attn_map = l.self_attn(x, mask=mask, return_attention=True)
            attention_maps.append(attn_map)
            x = l(x)
        return attention_maps
import math
from torch import Tensor

# from Moformer
# class PositionalEncoding(nn.Module):

#     def __init__(self, d_model: int, dropout: float = 0.1, max_len: int = 2048):
#         super().__init__()
#         self.dropout = nn.Dropout(p=dropout)

#         position = torch.arange(max_len).unsqueeze(1)
#         div_term = torch.exp(torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model))
#         pe = torch.zeros(max_len, 1, d_model)
#         pe[:, 0, 0::2] = torch.sin(position * div_term)
#         pe[:, 0, 1::2] = torch.cos(position * div_term)
#         self.register_buffer('pe', pe)

#     def forward(self, x: Tensor) -> Tensor:
#         """
#         Args:
#             x: Tensor, shape [seq_len, batch_size, embedding_dim]
#         """
#         x = x + self.pe[:x.size(0)]
#         return self.dropout(x)

# from UVDALC 
class PositionalEncoding(nn.Module):

    def __init__(self, d_model, max_len=5000):
        """
        Inputs
            d_model - Hidden dimensionality of the input.
            max_len - Maximum length of a sequence to expect.
        """
        super().__init__()

        # Create matrix of [SeqLen, HiddenDim] representing the positional encoding for max_len inputs
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)

        # register_buffer => Tensor which is not a parameter, but should be part of the modules state.
        # Used for tensors that need to be on the same device as the module.
        # persistent=False tells PyTorch to not add the buffer to the state dict (e.g. when we save the model)
        self.register_buffer('pe', pe, persistent=False)

    def forward(self, x):
        x = x + self.pe[:, :x.size(1)]
        return x
    
class TransformerModel(nn.Module):
    
    def __init__(self, args):
        super().__init__()
        num_layers=6
        input_dim=512
        model_dim=512
        dim_feedforward=2*512
        num_heads=8
        dropout=0
        num_classes=1
        input_dropout=0
        vocab_size = args['vocab_size']
        self.masked_token_train=args['masked_token_train']
        self.finetune = args['finetune']

#         with open(dictionary_path, 'rb') as f:
#             self.dictionary = pickle.load(f)
        
        self.embed = nn.Embedding(vocab_size, 512, 0)
        
        # self.input_net = nn.Sequential(
        #     nn.Dropout(input_dropout),
        #     nn.Linear(input_dim, hparams.model_dim)
        # )
        
        self.transformer = TransformerEncoder(num_layers=num_layers,
                                      input_dim=model_dim,
                                      dim_feedforward=2*model_dim,
                                      num_heads=num_heads,
                                      dropout=dropout)
        
        self.masked_model = MaskLMHead(embed_dim=512,
                                       output_dim=vocab_size
                                      )
        
        
        
        self.regression_head = ClassificationHead(input_dim=512,
                                                  inner_dim=256,
                                                  num_classes=1)
        
        self.output_net = nn.Sequential(
            nn.Linear(model_dim, model_dim),
            nn.LayerNorm(model_dim),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(model_dim, num_classes)
        )
            
     
        self.pos_encoder = PositionalEncoding(d_model=512, max_len=50)
            
    def forward(self, data ):
        attn_bias=None
        tokens = data.tokens
        
            
        """
        Inputs:
            x - Input features of shape [Batch, SeqLen, input_dim]
            mask - Mask to apply on the attention outputs (optional)
            add_positional_encoding - If True, we add the positional encoding to the input.
                                      Might not be desired for some tasks.
        """
        
        # tokens = batch['tokens']
        
        # masked_token_train=True
        token_mask=None
        if self.masked_token_train:
            
            rand = torch.rand(tokens.shape)
            # where the random array is less than 0.15, we set true
            mask_arr = rand < 0.15
            token_mask = mask_arr*(tokens !=1)*(tokens!=0)*(tokens!=2)

            tokens = tokens.masked_fill(token_mask.to(torch.bool), 4)

            

        
        
        padding_mask = tokens.eq(0)
        x = self.embed(tokens)
        x = self.pos_encoder(x) # adding positional encoding from Moformer paper https://github.com/zcao0420/MOFormer/blob/main/model/transformer.py
        # x = self.input_net(x)
        # if add_positional_encoding:
        #     x = self.positional_encoding(x)
        
        # encoder output
        encoder_out = self.transformer(x, mask=padding_mask, attn_bias=None)
        
        if self.masked_token_train:
            logits = self.masked_model(encoder_out, token_mask)
        else:
            logits = encoder_out
        
        
        if self.finetune:
            logits = self.regression_head(logits)
            
        
        
        # x = self.output_net(x)
        return logits, token_mask, encoder_out[:,0]
    
    
class MaskLMHead(nn.Module):
    """Head for masked language modeling."""

    def __init__(self, embed_dim, output_dim):
        super().__init__()
        self.dense = nn.Linear(embed_dim, embed_dim)
        self.activation_fn = F.gelu
        # self.layer_norm = LayerNorm(embed_dim)
        self.out = nn.Linear(embed_dim, output_dim)

        # if weight is None:
        #     weight = nn.Linear(embed_dim, output_dim, bias=False).weight
        # self.weight = weight
        # self.bias = nn.Parameter(torch.zeros(output_dim))

    def forward(self, features, masked_tokens=None, **kwargs):
        # Only project the masked tokens while training,
        # saves both memory and computation
        if masked_tokens is not None:
            features = features[masked_tokens, :]

        x = self.dense(features)
        x = self.activation_fn(x)
        # x = self.layer_norm(x)
        # project back to size of vocabulary with bias
        x = self.out(x)
        return x


class ClassificationHead(nn.Module):
    """Head for sentence-level classification tasks."""

    def __init__(
        self,
        input_dim,
        inner_dim,
        num_classes,
        pooler_dropout=0.0,
    ):
        super().__init__()
        self.dense = nn.Linear(input_dim, inner_dim)
        self.activation_fn = F.relu
        self.dropout = nn.Dropout(p=pooler_dropout)
        self.out_proj = nn.Linear(inner_dim, num_classes)

    def forward(self, features, **kwargs):
        x = features[:, 0, :]  # take <s> token (equiv. to [CLS])
        x = self.dropout(x)
        x = self.dense(x)
        x = self.activation_fn(x)
        x = self.dropout(x)
        x = self.out_proj(x)
        return x
