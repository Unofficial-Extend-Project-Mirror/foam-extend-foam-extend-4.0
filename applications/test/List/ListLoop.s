	.file	"ListLoop.C"
	.section	.ctors,"aw",@progbits
	.align 4
	.long	_GLOBAL__I__Z4funcRN4Foam4ListIdEERKS1_RKNS0_IiEES7_
	.text
	.align 2
	.p2align 4,,15
	.type	_Z41__static_initialization_and_destruction_0ii, @function
_Z41__static_initialization_and_destruction_0ii:
.LFB2550:
	pushl	%ebp
.LCFI0:
	movl	%esp, %ebp
.LCFI1:
	pushl	%ebx
.LCFI2:
	call	__i686.get_pc_thunk.bx
	addl	$_GLOBAL_OFFSET_TABLE_, %ebx
	subl	$20, %esp
.LCFI3:
	decl	%eax
	je	.L7
.L5:
	addl	$20, %esp
	popl	%ebx
	leave
	ret
	.p2align 4,,7
.L7:
	cmpl	$65535, %edx
	jne	.L5
	leal	_ZSt8__ioinit@GOTOFF(%ebx), %eax
	movl	%eax, (%esp)
	call	_ZNSt8ios_base4InitC1Ev@PLT
	movl	__dso_handle@GOT(%ebx), %eax
	movl	$0, 4(%esp)
	movl	%eax, 8(%esp)
	leal	__tcf_0@GOTOFF(%ebx), %eax
	movl	%eax, (%esp)
	call	__cxa_atexit@PLT
	addl	$20, %esp
	popl	%ebx
	leave
	ret
.LFE2550:
	.size	_Z41__static_initialization_and_destruction_0ii, .-_Z41__static_initialization_and_destruction_0ii
.globl __gxx_personality_v0
	.align 2
	.p2align 4,,15
	.type	_GLOBAL__I__Z4funcRN4Foam4ListIdEERKS1_RKNS0_IiEES7_, @function
_GLOBAL__I__Z4funcRN4Foam4ListIdEERKS1_RKNS0_IiEES7_:
.LFB2552:
	pushl	%ebp
.LCFI4:
	movl	$65535, %edx
	movl	$1, %eax
	movl	%esp, %ebp
.LCFI5:
	leave
	jmp	_Z41__static_initialization_and_destruction_0ii
.LFE2552:
	.size	_GLOBAL__I__Z4funcRN4Foam4ListIdEERKS1_RKNS0_IiEES7_, .-_GLOBAL__I__Z4funcRN4Foam4ListIdEERKS1_RKNS0_IiEES7_
	.align 2
	.p2align 4,,15
	.type	__tcf_0, @function
__tcf_0:
.LFB2551:
	pushl	%ebp
.LCFI6:
	movl	%esp, %ebp
.LCFI7:
	pushl	%ebx
.LCFI8:
	call	__i686.get_pc_thunk.bx
	addl	$_GLOBAL_OFFSET_TABLE_, %ebx
	subl	$4, %esp
.LCFI9:
	leal	_ZSt8__ioinit@GOTOFF(%ebx), %eax
	movl	%eax, (%esp)
	call	_ZNSt8ios_base4InitD1Ev@PLT
	addl	$4, %esp
	popl	%ebx
	leave
	ret
.LFE2551:
	.size	__tcf_0, .-__tcf_0
	.align 2
	.p2align 4,,15
.globl _Z4funcRN4Foam4ListIdEERKS1_RKNS0_IiEES7_
	.type	_Z4funcRN4Foam4ListIdEERKS1_RKNS0_IiEES7_, @function
_Z4funcRN4Foam4ListIdEERKS1_RKNS0_IiEES7_:
.LFB2352:
	pushl	%ebp
.LCFI10:
	movl	%esp, %ebp
.LCFI11:
	pushl	%edi
.LCFI12:
	pushl	%esi
.LCFI13:
	subl	$12, %esp
.LCFI14:
	movl	8(%ebp), %edx
	movl	(%edx), %eax
	testl	%eax, %eax
	movl	%eax, -20(%ebp)
	jle	.L16
	movl	16(%ebp), %eax
	movl	4(%edx), %edx
	xorl	%ecx, %ecx
	movl	4(%eax), %eax
	movl	%edx, -12(%ebp)
	movl	%eax, -16(%ebp)
	movl	20(%ebp), %eax
	movl	4(%eax), %edi
	movl	12(%ebp), %eax
	movl	4(%eax), %esi
	.p2align 4,,7
.L15:
	movl	-16(%ebp), %edx
	movl	(%edx,%ecx,4), %eax
	movl	-12(%ebp), %edx
	leal	(%edx,%eax,8), %eax
	movl	(%edi,%ecx,4), %edx
	incl	%ecx
	cmpl	%ecx, -20(%ebp)
	fldl	(%eax)
	fsubl	(%esi,%edx,8)
	fstpl	(%eax)
	jne	.L15
.L16:
	addl	$12, %esp
	popl	%esi
	popl	%edi
	leave
	ret
.LFE2352:
	.size	_Z4funcRN4Foam4ListIdEERKS1_RKNS0_IiEES7_, .-_Z4funcRN4Foam4ListIdEERKS1_RKNS0_IiEES7_
	.local	_ZSt8__ioinit
	.comm	_ZSt8__ioinit,1,1
	.section	.eh_frame,"a",@progbits
.Lframe1:
	.long	.LECIE1-.LSCIE1
.LSCIE1:
	.long	0x0
	.byte	0x1
	.string	"zPR"
	.uleb128 0x1
	.sleb128 -4
	.byte	0x8
	.uleb128 0x6
	.byte	0x9b
	.long	DW.ref.__gxx_personality_v0-.
	.byte	0x1b
	.byte	0xc
	.uleb128 0x4
	.uleb128 0x4
	.byte	0x88
	.uleb128 0x1
	.align 4
.LECIE1:
.LSFDE1:
	.long	.LEFDE1-.LASFDE1
.LASFDE1:
	.long	.LASFDE1-.Lframe1
	.long	.LFB2550-.
	.long	.LFE2550-.LFB2550
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI0-.LFB2550
	.byte	0xe
	.uleb128 0x8
	.byte	0x85
	.uleb128 0x2
	.byte	0x4
	.long	.LCFI1-.LCFI0
	.byte	0xd
	.uleb128 0x5
	.byte	0x4
	.long	.LCFI2-.LCFI1
	.byte	0x83
	.uleb128 0x3
	.align 4
.LEFDE1:
.LSFDE5:
	.long	.LEFDE5-.LASFDE5
.LASFDE5:
	.long	.LASFDE5-.Lframe1
	.long	.LFB2551-.
	.long	.LFE2551-.LFB2551
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI6-.LFB2551
	.byte	0xe
	.uleb128 0x8
	.byte	0x85
	.uleb128 0x2
	.byte	0x4
	.long	.LCFI7-.LCFI6
	.byte	0xd
	.uleb128 0x5
	.byte	0x4
	.long	.LCFI8-.LCFI7
	.byte	0x83
	.uleb128 0x3
	.align 4
.LEFDE5:
	.hidden DW.ref.__gxx_personality_v0
	.weak	DW.ref.__gxx_personality_v0
	.section	.data.DW.ref.__gxx_personality_v0,"awG",@progbits,DW.ref.__gxx_personality_v0,comdat
	.align 4
	.type	DW.ref.__gxx_personality_v0, @object
	.size	DW.ref.__gxx_personality_v0, 4
DW.ref.__gxx_personality_v0:
	.long	__gxx_personality_v0
	.ident	"GCC: (GNU) 4.1.1"
	.section	.text.__i686.get_pc_thunk.bx,"axG",@progbits,__i686.get_pc_thunk.bx,comdat
.globl __i686.get_pc_thunk.bx
	.hidden	__i686.get_pc_thunk.bx
	.type	__i686.get_pc_thunk.bx, @function
__i686.get_pc_thunk.bx:
	movl	(%esp), %ebx
	ret
	.section	.note.GNU-stack,"",@progbits
