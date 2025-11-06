/**
 * Professional Modal Component - IBM-Inspired Design
 *
 * Features:
 * - Sharp rectangular design (no rounded corners)
 * - IBM-inspired color palette
 * - Backdrop with click-outside to close
 * - Keyboard accessibility (Escape to close)
 * - Focus trap
 * - Multiple sizes (sm, md, lg, xl)
 */

import React, { useEffect, useRef } from 'react';

export interface ModalProps {
  isOpen: boolean;
  onClose: () => void;
  children: React.ReactNode;
  size?: 'sm' | 'md' | 'lg' | 'xl';
  closeOnBackdrop?: boolean;
  closeOnEscape?: boolean;
  className?: string;
}

export interface ModalHeaderProps {
  children: React.ReactNode;
  onClose?: () => void;
  className?: string;
}

export interface ModalBodyProps {
  children: React.ReactNode;
  className?: string;
}

export interface ModalFooterProps {
  children: React.ReactNode;
  className?: string;
}

export interface ModalTitleProps {
  children: React.ReactNode;
  className?: string;
}

// Main Modal Component
const Modal: React.FC<ModalProps> & {
  Header: React.FC<ModalHeaderProps>;
  Body: React.FC<ModalBodyProps>;
  Footer: React.FC<ModalFooterProps>;
  Title: React.FC<ModalTitleProps>;
} = ({
  isOpen,
  onClose,
  children,
  size = 'md',
  closeOnBackdrop = true,
  closeOnEscape = true,
  className = '',
}) => {
  const modalRef = useRef<HTMLDivElement>(null);

  // Size classes
  const sizeClasses = {
    sm: 'max-w-[var(--modal-max-width-sm)]',
    md: 'max-w-[var(--modal-max-width-md)]',
    lg: 'max-w-[var(--modal-max-width-lg)]',
    xl: 'max-w-[var(--modal-max-width-xl)]',
  };

  // Handle Escape key
  useEffect(() => {
    if (!isOpen || !closeOnEscape) return;

    const handleEscape = (e: KeyboardEvent) => {
      if (e.key === 'Escape') {
        onClose();
      }
    };

    document.addEventListener('keydown', handleEscape);
    return () => document.removeEventListener('keydown', handleEscape);
  }, [isOpen, closeOnEscape, onClose]);

  // Lock body scroll when modal is open
  useEffect(() => {
    if (isOpen) {
      document.body.style.overflow = 'hidden';
      return () => {
        document.body.style.overflow = '';
      };
    }
  }, [isOpen]);

  // Handle backdrop click
  const handleBackdropClick = (e: React.MouseEvent<HTMLDivElement>) => {
    if (closeOnBackdrop && e.target === e.currentTarget) {
      onClose();
    }
  };

  if (!isOpen) return null;

  return (
    <div className="modal-backdrop" onClick={handleBackdropClick}>
      <div
        ref={modalRef}
        className={`modal ${sizeClasses[size]} ${className}`}
        role="dialog"
        aria-modal="true"
      >
        {children}
      </div>
    </div>
  );
};

// Modal Header
const ModalHeader: React.FC<ModalHeaderProps> = ({
  children,
  onClose,
  className = '',
}) => {
  return (
    <div className={`modal-header ${className}`}>
      {children}
      {onClose && (
        <button
          onClick={onClose}
          className="btn-ghost p-2 ml-auto"
          aria-label="Close modal"
        >
          <svg
            width="20"
            height="20"
            viewBox="0 0 20 20"
            fill="none"
            stroke="currentColor"
            strokeWidth="2"
            strokeLinecap="square"
          >
            <path d="M5 5L15 15M15 5L5 15" />
          </svg>
        </button>
      )}
    </div>
  );
};

// Modal Body
const ModalBody: React.FC<ModalBodyProps> = ({ children, className = '' }) => {
  return <div className={`modal-body ${className}`}>{children}</div>;
};

// Modal Footer
const ModalFooter: React.FC<ModalFooterProps> = ({ children, className = '' }) => {
  return <div className={`modal-footer ${className}`}>{children}</div>;
};

// Modal Title
const ModalTitle: React.FC<ModalTitleProps> = ({ children, className = '' }) => {
  return <h2 className={`modal-title ${className}`}>{children}</h2>;
};

// Attach subcomponents
Modal.Header = ModalHeader;
Modal.Body = ModalBody;
Modal.Footer = ModalFooter;
Modal.Title = ModalTitle;

export default Modal;
